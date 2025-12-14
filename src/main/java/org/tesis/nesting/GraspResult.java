package org.tesis.nesting;

import java.util.List;

/**
 * Clase para almacenar el resultado final del algoritmo GRASP,
 * incluyendo la solución y las estadísticas.
 */
public class GraspResult {
    public List<PlacedPart> solution;
    public int corrida;
    public long seed;
    public int toPlace;
    public int placedCount;
    public double areaPlaced;
    public double areaSheet;
    public double areaWaste;
    public double yMax;
    public double utilizationPct;
    public long elapsedMs;
    public String algo = "GRASP";
}


