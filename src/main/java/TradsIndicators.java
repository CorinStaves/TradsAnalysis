import com.google.common.math.LongMath;
import org.apache.log4j.Logger;
import org.geotools.data.simple.SimpleFeatureReader;
import org.geotools.geopkg.GeoPackage;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;
import org.matsim.api.core.v01.Coord;
import org.matsim.api.core.v01.TransportMode;
import org.matsim.api.core.v01.network.*;
import org.matsim.api.core.v01.population.Person;
import org.matsim.core.config.Config;
import org.matsim.core.config.ConfigUtils;
import org.matsim.core.network.NetworkUtils;
import org.matsim.core.network.algorithms.TransportModeNetworkFilter;
import org.matsim.core.network.io.MatsimNetworkReader;
import org.matsim.core.router.FastDijkstraFactory;
import org.matsim.core.router.costcalculators.FreespeedTravelTimeAndDisutility;
import org.matsim.core.router.util.LeastCostPathCalculator;
import org.matsim.core.router.util.TravelDisutility;
import org.matsim.core.router.util.TravelTime;
import org.matsim.core.utils.geometry.CoordUtils;
import org.matsim.core.utils.misc.Counter;
import org.matsim.vehicles.Vehicle;
import org.opengis.feature.simple.SimpleFeature;

import java.io.*;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.function.Predicate;

// THIS SCRIPT IS FOR CONVERTING X/Y COORDINATES IN THE TRAVEL SURVEY TO PATH DATA (E.G. TRAVEL TIMES, DISTANCES, COSTS) USING THE MATSIM NETWORK
// Currently missing public transport

public class TradsIndicators {

    private final static Logger logger = Logger.getLogger(TradsIndicators.class);

    private final static String HOUSEHOLD_ID = "IDNumber";
    private final static String PERSON_ID = "PersonNumber";
    private final static String TRIP_ID = "TripNumber";

    private final static String SEP = ",";

    private final static String X_ORIGIN_COORD = "StartEasting";
    private final static String Y_ORIGIN_COORD = "StartNorthing";

    private final static String X_DESTINATION_COORD = "EndEasting";
    private final static String Y_DESTINATION_COORD = "EndNorthing";

    private final static double WALK_SPEED = 5.3 / 3.6;
    private final static double BIKE_SPEED = 14 / 3.6;

    public static void main(String[] args) throws IOException {

        if(args.length != 5) {
            throw new RuntimeException("Program requires 5 arguments: \n" +
                    "(1) Survey File Path \n" +
                    "(2) Boundary Geopackage Path \n" +
                    "(3) Network File Path \n" +
                    "(4) Output File Path \n" +
                    "(5) Number of Threads \n");
        }

        String surveyFilePath = args[0];
        String boundaryFilePath = args[1];
        String networkFilePath = args[2];
        String outputFile = args[3];
        int numberOfThreads = Integer.parseInt(args[4]);

        // Read network
        logger.info("Reading MATSim network...");
        Network network = NetworkUtils.createNetwork();
        new MatsimNetworkReader(network).readFile(networkFilePath);

        // Set up scenario and config
        logger.info("Preparing Matsim config and scenario...");
        Config config = ConfigUtils.createConfig();

        // Create mode-specific networks
        logger.info("Creating mode-specific networks...");
        Network networkCar = extractModeSpecificNetwork(network, TransportMode.car);
        Network carXy2l = extractXy2LinksNetwork(networkCar, l -> !((boolean) l.getAttributes().getAttribute("motorway")));
        Network networkBike = extractModeSpecificNetwork(network, TransportMode.bike);
        Network networkWalk = extractModeSpecificNetwork(network, TransportMode.walk);

        // Read Boundary Shapefile
        logger.info("Reading boundary shapefile...");
        Geometry boundary = readNetworkBoundary(boundaryFilePath);

        // Read in coordinates from CSV
        logger.info("Reading person micro data from ascii file...");
        Set<Trip> trips = readCoordinates(surveyFilePath, boundary);

        // Travel time and disutility data
        FreespeedTravelTimeAndDisutility freeSpeed = new FreespeedTravelTimeAndDisutility(config.planCalcScore());
        TravelDisutility distanceAsTravelDisutility = new DistanceAsTravelDisutility();
        TravelTime travelTimeBike = (link, v, person, vehicle) -> link.getLength() / BIKE_SPEED;
        TravelTime travelTimeWalk = (link, v, person, vehicle) -> link.getLength() / WALK_SPEED;

        // Calculate indicators
        logger.info("Calculating network indicators using " + numberOfThreads + " threads.");
        calculateNetworkIndicators(trips, "car", networkCar, carXy2l, numberOfThreads, freeSpeed, freeSpeed);
        calculateNetworkIndicators(trips, "bike", networkBike, null, numberOfThreads, distanceAsTravelDisutility, travelTimeBike);
        calculateNetworkIndicators(trips, "walk", networkWalk, null, numberOfThreads, distanceAsTravelDisutility, travelTimeWalk);

        // Write results
        logger.info("Writing results to csv file...");
        writeIndicators(trips, outputFile);
    }

    private static Geometry readNetworkBoundary(String filePath) throws IOException {
        GeoPackage geopkg = new GeoPackage(new File(filePath));
        SimpleFeatureReader r = geopkg.reader(geopkg.features().get(0), null,null);
        SimpleFeature f = r.next();
        Geometry boundary = (Geometry) f.getDefaultGeometry();
        r.close();
        geopkg.close();
        return boundary;
    }

    private static int findPositionInArray (String string, String[] array) {
        int ind = -1;
        for (int a = 0; a < array.length; a++) {
            if (array[a].equalsIgnoreCase(string)) {
                ind = a;
            }
        }
        if (ind == -1) {
            logger.error ("Could not find element " + string +
                    " in array (see method <findPositionInArray> in class <SiloUtil>");
        }
        return ind;
    }

    private static PrintWriter openFileForSequentialWriting(File outputFile) {
        if (outputFile.getParent() != null) {
            File parent = outputFile.getParentFile();
            parent.mkdirs();
        }

        try {
            FileWriter fw = new FileWriter(outputFile, false);
            BufferedWriter bw = new BufferedWriter(fw);
            return new PrintWriter(bw);
        } catch (IOException var5) {
            logger.info("Could not open file <" + outputFile.getName() + ">.");
            return null;
        }
    }

    // Extracts mode-specific network  (e.g. walk network, car network, cycle network)
    private static Network extractModeSpecificNetwork(Network network, String transportMode) {
        Network modeSpecificNetwork = NetworkUtils.createNetwork();
        new TransportModeNetworkFilter(network).filter(modeSpecificNetwork, Collections.singleton(transportMode));
        NetworkUtils.runNetworkCleaner(modeSpecificNetwork);
        return modeSpecificNetwork;
    }

    // Extracts network of usable nearest links to start/end journey (e.g. a car trip cannot start on a motorway)
    private static Network extractXy2LinksNetwork(Network network, Predicate<Link> xy2linksPredicate) {
        Network xy2lNetwork = NetworkUtils.createNetwork();
        NetworkFactory nf = xy2lNetwork.getFactory();
        for (Link link : network.getLinks().values()) {
            if (xy2linksPredicate.test(link)) {
                // okay, we need that link
                Node fromNode = link.getFromNode();
                Node xy2lFromNode = xy2lNetwork.getNodes().get(fromNode.getId());
                if (xy2lFromNode == null) {
                    xy2lFromNode = nf.createNode(fromNode.getId(), fromNode.getCoord());
                    xy2lNetwork.addNode(xy2lFromNode);
                }
                Node toNode = link.getToNode();
                Node xy2lToNode = xy2lNetwork.getNodes().get(toNode.getId());
                if (xy2lToNode == null) {
                    xy2lToNode = nf.createNode(toNode.getId(), toNode.getCoord());
                    xy2lNetwork.addNode(xy2lToNode);
                }
                Link xy2lLink = nf.createLink(link.getId(), xy2lFromNode, xy2lToNode);
                xy2lLink.setAllowedModes(link.getAllowedModes());
                xy2lLink.setCapacity(link.getCapacity());
                xy2lLink.setFreespeed(link.getFreespeed());
                xy2lLink.setLength(link.getLength());
                xy2lLink.setNumberOfLanes(link.getNumberOfLanes());
                xy2lNetwork.addLink(xy2lLink);
            }
        }
        return xy2lNetwork;
    }

    // Reads Survey File with Coordinate Data
    private static Set<Trip> readCoordinates(String filePath, Geometry geometry) throws IOException {
        Set<Trip> trips = new HashSet<>();
        String recString = "";
        int counter = 0;
        int badRecords = 0;
        int originsOutsideBoundary = 0;
        int destinationsOutsideBoundary = 0;

        // Open Reader
        BufferedReader in = new BufferedReader(new FileReader(filePath));

        GeometryFactory gf = new GeometryFactory();

        // Read Header
        recString = in.readLine();
        String[] header = recString.split(SEP);
        int posHouseholdId = findPositionInArray(HOUSEHOLD_ID, header);
        int posPersonId = findPositionInArray(PERSON_ID, header);
        int posTripId = findPositionInArray(TRIP_ID, header);
        int posOrigX = findPositionInArray(X_ORIGIN_COORD, header);
        int posOrigY = findPositionInArray(Y_ORIGIN_COORD, header);
        int posDestX = findPositionInArray(X_DESTINATION_COORD, header);
        int posDestY = findPositionInArray(Y_DESTINATION_COORD, header);

        while ((recString = in.readLine()) != null) {
            counter++;
            if (LongMath.isPowerOfTwo(counter)) {
                logger.info(counter + " records processed.");
            }
            String[] lineElements = recString.split(SEP);

            String householdId = lineElements[posHouseholdId];
            int personId = Integer.parseInt(lineElements[posPersonId]);
            int tripId = Integer.parseInt(lineElements[posTripId]);
            Coord cOrig = null;
            Coord cDest = null;
            Boolean originInBoundary = null;
            Boolean destinationInBoundary = null;
            try {
                double xOrig = Double.parseDouble(lineElements[posOrigX]);
                double yOrig = Double.parseDouble(lineElements[posOrigY]);
                double xDest = Double.parseDouble(lineElements[posDestX]);
                double yDest = Double.parseDouble(lineElements[posDestY]);

                cOrig = CoordUtils.createCoord(xOrig, yOrig);
                cDest = CoordUtils.createCoord(xDest, yDest);

                originInBoundary = geometry.contains(gf.createPoint(new Coordinate(xOrig, yOrig)));
                destinationInBoundary = geometry.contains(gf.createPoint(new Coordinate(xDest, yDest)));

                if (!originInBoundary) originsOutsideBoundary++;
                if (!destinationInBoundary) destinationsOutsideBoundary++;

            } catch (NumberFormatException e) {
                logger.warn("Unreadable coordinates for record at line " + counter);
                badRecords++;
            }

            trips.add(new Trip(householdId, personId, tripId, cOrig, cDest, originInBoundary, destinationInBoundary));
        }
        in.close();
        logger.info(counter + " records processed.");
        logger.info(trips.size() + " trips read successfuly.");
        logger.info(originsOutsideBoundary + " trips with origin coordinates outside study area boundary.");
        logger.info(destinationsOutsideBoundary + " trips with destination coordinates outside study area boundary.");
        logger.info(badRecords + " trips with no or unreadable coordinates.");

        return trips;
    }

    // Writes results to file
    private static void writeIndicators(Set<Trip> trips, String filePath) {

        PrintWriter out = openFileForSequentialWriting(new File(filePath));
        out.println(HOUSEHOLD_ID + SEP + PERSON_ID + SEP + TRIP_ID + SEP +
                "OriginWithinBoundary" + SEP + "DestinationWithinBoundary" + SEP +
                "CarCost" + SEP + "CarTime" + SEP +
                "BikeCost" + SEP + "BikeTime" + SEP +
                "WalkCost" + SEP + "WalkTime");

        int counter = 0;
        for (Trip trip : trips) {
            counter++;
            if (LongMath.isPowerOfTwo(counter)) {
                logger.info(counter + " records written.");
            }
            out.println(trip.getHouseholdId() + SEP + trip.getPersonId() + SEP + trip.getTripId() + SEP +
                    trip.isOriginWithinBoundary() + SEP + trip.isDestinationWithinBoundary() + SEP +
                    trip.getCost("car") + SEP + trip.getTime("car") + SEP +
                    trip.getCost("bike") + SEP + trip.getTime("bike") + SEP +
                    trip.getCost("walk") + SEP + trip.getTime("walk"));
        }
        out.close();
        logger.info("Wrote " + counter + " records to " + filePath);
    }

    private static class DistanceAsTravelDisutility implements TravelDisutility {
        public DistanceAsTravelDisutility() {
        }

        public double getLinkTravelDisutility(Link link, double time, Person person, Vehicle vehicle) {
            return link.getLength();
        }

        public double getLinkMinimumTravelDisutility(Link link) {
            return link.getLength();
        }
    }

    // Calculates cost (km) and time (seconds) for each mode & network (multithreaded)
    private static void calculateNetworkIndicators(Set<Trip> trips, String mode, Network network, Network xy2lNetwork, int numberOfThreads,
                                                   TravelDisutility travelDisutility, TravelTime travelTime) {

        logger.info("Calculating indicators for mode " + mode);

        // do calculation
        ConcurrentLinkedQueue<Trip> odPairsQueue = new ConcurrentLinkedQueue<>(trips);

        Counter counter = new Counter(mode + ": Route ", " / " + trips.size());
        Thread[] threads = new Thread[numberOfThreads];
        for (int i = 0; i < numberOfThreads; i++) {
            NetworkIndicatorCalculator worker = new NetworkIndicatorCalculator(odPairsQueue, counter, mode, network, xy2lNetwork, travelDisutility, travelTime);
            threads[i] = new Thread(worker, "IndicatorCalculator-" + mode + "-" + i);
            threads[i].start();
        }

        // wait until all threads have finished
        for (Thread thread : threads) {
            try {
                thread.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    private static class NetworkIndicatorCalculator implements Runnable {

        private final ConcurrentLinkedQueue<Trip> trips;
        private final Counter counter;
        private final String mode;
        private final Network routingNetwork;
        private final Network xy2lNetwork;
        private final TravelDisutility travelDisutility;
        private final TravelTime travelTime;

        NetworkIndicatorCalculator(ConcurrentLinkedQueue<Trip> trips, Counter counter, String mode,
                                   Network routingNetwork, Network xy2lNetwork,
                                   TravelDisutility travelDisutility, TravelTime travelTime) {
            this.trips = trips;
            this.counter = counter;
            this.mode = mode;
            this.routingNetwork = routingNetwork;
            this.xy2lNetwork = xy2lNetwork;
            this.travelDisutility = travelDisutility;
            this.travelTime = travelTime;
        }

        public void run() {
            LeastCostPathCalculator dijkstra = new FastDijkstraFactory(false).
                    createPathCalculator(this.routingNetwork, this.travelDisutility, this.travelTime);

            while(true) {
                Trip trip = this.trips.poll();
                if(trip == null) {
                    return;
                }

                this.counter.incCounter();

                if(trip.isTripWithinBoundary()) {
                    Coord cOrig = trip.getOrigCoord();
                    Coord cDest = trip.getDestCoord();
                    Node nOrig;
                    Node nDest;
                    if(xy2lNetwork == null) {
                        nOrig = NetworkUtils.getNearestNode(routingNetwork,cOrig);
                        nDest = NetworkUtils.getNearestNode(routingNetwork,cDest);
                    } else {
                        nOrig = routingNetwork.getNodes().get(NetworkUtils.getNearestLink(xy2lNetwork, cOrig).getToNode().getId());
                        nDest = routingNetwork.getNodes().get(NetworkUtils.getNearestLink(xy2lNetwork, cDest).getToNode().getId());
                    }

                    // Calculate least cost path
                    LeastCostPathCalculator.Path path = dijkstra.calcLeastCostPath(nOrig, nDest, 28800, null, null);
                    double cost = path.travelCost / 1000;
                    double time = path.travelTime;

                    // Store processed data
                    trip.addProcessedData(mode, cost, time);
                } else {
                    // Cannot calculate if trip is not within boundary
                    trip.addProcessedData(mode, Double.NaN, Double.NaN);
                }
            }
        }
    }
}
