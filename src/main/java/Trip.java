import org.matsim.api.core.v01.Coord;

import java.util.HashMap;
import java.util.Map;

public class Trip {

    private final String householdId;
    private final int personId;
    private final int tripId;
    
    private final Coord cOrig;
    private final Coord cDest;

    private Map<String, Double> costs = new HashMap<>();
    private Map<String, Double> travelTimes = new HashMap<>();

    private Boolean originWitinBoundary;
    private Boolean destinationWithinBoundary;

    public Trip(String householdId, int personId, int tripId, Coord cOrig, Coord cDest,
                Boolean originWithinBoundary, Boolean destinationWithinBoundary) {
        this.householdId = householdId;
        this.personId = personId;
        this.tripId = tripId;
        this.cOrig = cOrig;
        this.cDest = cDest;
        this.originWitinBoundary = originWithinBoundary;
        this.destinationWithinBoundary = destinationWithinBoundary;
    }

    public Boolean isOriginWithinBoundary() {
        return originWitinBoundary;
    }

    public Boolean isDestinationWithinBoundary() {
        return destinationWithinBoundary;
    }

    public boolean isTripWithinBoundary() {
        if(originWitinBoundary != null && destinationWithinBoundary != null) {
            return originWitinBoundary && destinationWithinBoundary;
        }
        else return false;
    }

    public Coord getOrigCoord() {
        return cOrig;
    }
    
    public Coord getDestCoord() {
        return cDest;
    }
    
    public void addProcessedData(String mode, double cost, double time) {
        costs.put(mode, cost);
        travelTimes.put(mode, time);
    }

    public String getHouseholdId() {
        return householdId;
    }

    public int getPersonId() {
        return personId;
    }

    public int getTripId() {
        return tripId;
    }

    public double getCost(String mode) {
        return costs.get(mode);
    }

    public double getTime(String mode) {
        return travelTimes.get(mode);
    }
    
}
