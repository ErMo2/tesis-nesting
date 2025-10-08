package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.operation.union.UnaryUnionOp;

import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * =====================  HybridGraspMemeticPlacer  =====================
 * Estrategia híbrida:
 *   1) Corre GRASP con múltiples semillas para generar soluciones iniciales.
 *   2) Pasa esas soluciones al memético (HybridMemetic) para refinarlas.
 *   3) Devuelve la mejor solución final.
 *
 */
public class HybridGraspMemeticPlacer {

    // ===================== Config Híbrido =====================
    private static final int    GRASP_SEEDS       = 5;          // # de corridas GRASP
    private static final long   MEMETIC_TIME_MS   = 120_000L;   // tiempo máx. memético (ms)
    private static final String LOG_DIR           = "";
    private static final String LOG_FILE          = File.separator + "LogHybrid.log";

    // Formatos
    private static final DecimalFormat DF3  = new DecimalFormat("#,##0.000");
    private static final SimpleDateFormat TS = new SimpleDateFormat("HH:mm:ss.SSS");

    public static List<PlacedPart> place(Polygon sheet, List<PartSpec> parts) throws Exception {
        return place(sheet, parts, GRASP_SEEDS, MEMETIC_TIME_MS, LOG_FILE);
    }

    // Overload con parámetros (por si quieres tunear desde tu main)
    public static List<PlacedPart> place(Polygon sheet,
                                         List<PartSpec> parts,
                                         int graspSeeds,
                                         long memeticTimeMs,
                                         String logFile) throws Exception {
        
        Objects.requireNonNull(sheet, "sheet no puede ser null");
        Objects.requireNonNull(parts, "parts no puede ser null");
        new File(LOG_DIR).mkdirs();
        final int TOTAL_PARTS = parts.stream()
        .mapToInt(p -> Math.max(1, p.qty))   // suma todos los qty, usa 1 si es 0
        .sum();
        long seed0 = System.nanoTime();
        long t0 = System.currentTimeMillis();

        try (PrintStream log = openLog(logFile)) {

            head(log, "==== HYB GRASP→MEMETIC | START ====");
            line(log, "CFG | GRASP_SEEDS=" + graspSeeds + "  MEMETIC_TIME_MS=" + memeticTimeMs);
            line(log, "CFG | seed0=" + seed0 + "  parts=" + parts.size());
            
            
            
            // ============= 1) Multi-GRASP: generar población semilla ============
            List<List<PlacedPart>> seedPopulation = new ArrayList<>();
            List<Stats> graspStatsList = new ArrayList<>();

            for (int i = 0; i < graspSeeds; i++) {
                long seed = seed0 + i;
                long tg0 = System.currentTimeMillis();
                line(log, "GRASP seed " + String.format("%02d", i) + " (seed=" + seed + ") | start");
                List<PlacedPart> sol = GraspPlacer.place(sheet, new ArrayList<>(parts), seed);
                long tg1 = System.currentTimeMillis();

                Stats s = evalSolution("GRASP seed " + String.format("%02d", i),
                       sheet, TOTAL_PARTS, sol, tg1 - tg0, seed);
                graspStatsList.add(s);
                printStats(log, s);

                if (sol != null && !sol.isEmpty()) {
                    // sanitizar para asegurar validez antes de inyectar
                    seedPopulation.add(sanitizeSolution(sheet, sol));
                }
            }

            // Si por algún motivo no hubo ninguna sol válida, retornamos la mejor de GRASP (o vacío)
            if (seedPopulation.isEmpty()) {
                line(log, "WARN | No hubo soluciones válidas de GRASP para inyectar al memético. Devolviendo la mejor de GRASP (si existe).");
                Stats bestG = bestOf(graspStatsList);
                summary(log, "RESUMEN FINAL (HYBRID - solo GRASP)");
                if (bestG != null) {
                    renderSummary(log, TOTAL_PARTS, bestG, System.currentTimeMillis() - t0);
                    dash(log);
                    line(log, "BEST | GRASP-only | colocadas=" + bestG.placedCount +
                               " | area=" + DF3.format(bestG.areaPlaced) +
                               " | yMax=" + DF3.format(bestG.yMax));
                    dash(log);
                    return bestG.solution != null ? bestG.solution : Collections.emptyList();
                } else {
                    dash(log);
                    line(log, "BEST | vacío (no se pudo colocar nada)");
                    dash(log);
                    return Collections.emptyList();
                }
            }

            // ============= 2) Memético: refinar población semilla ===============
            line(log, "MEMETIC | improveFromPopulation | start");
            long tm0 = System.currentTimeMillis();
            List<PlacedPart> improved = HybridMemetic.improveFromPopulation(
                    sheet,
                    parts,
                    seedPopulation,
                    logFile,           // mismo log del híbrido
                    memeticTimeMs,
                    seed0 + 777       // semilla distinta para el memético
            );
            long tm1 = System.currentTimeMillis();

            Stats memStats = evalSolution("HYB-MEM", sheet, TOTAL_PARTS, improved, tm1 - tm0, seed0 + 777);
            printStats(log, memStats);

            // ============= 3) Seleccionar resultado final =======================
            // Compara mejor de GRASP vs mejor Memético (por si GRASP fue mejor)
            Stats bestGrasp = bestOf(graspStatsList);
            Stats finalStats = (compare(memStats, bestGrasp) < 0) ? memStats : bestGrasp;
            List<PlacedPart> finalSol = finalStats.solution;

            // ============= 4) Resumen final en el formato solicitado ============
            long totalMs = System.currentTimeMillis() - t0;
            summary(log, "RESUMEN FINAL (HYBRID)");
            renderSummary(log, parts.size(), finalStats, totalMs);
            dash(log);
            line(log, "BEST | " + finalStats.algo +
                       " | colocadas=" + finalStats.placedCount +
                       " | area=" + DF3.format(finalStats.areaPlaced) +
                       " | yMax=" + DF3.format(finalStats.yMax));
            dash(log);

            return (finalSol != null) ? finalSol : Collections.emptyList();
        }
    }

    // ===================== Métricas / evaluación =====================

    private static Stats evalSolution(String algoName,
                                  Polygon sheet,
                                  int totalParts,
                                  List<PlacedPart> placed,
                                  long elapsedMs,
                                  long seed) {
        Stats s = new Stats();
        s.algo = algoName;
        s.seed = seed;
        s.solution = placed;
        s.toPlace = totalParts; // ✅ ahora recibe el valor correcto
        s.placedCount = (placed != null) ? placed.size() : 0;
        s.areaSheet = sheet.getArea();
        s.areaPlaced = placedArea(placed);
        s.areaWaste = Math.max(0.0, s.areaSheet - s.areaPlaced);
        s.utilizationPct = (s.areaSheet <= 0) ? 0.0 : (s.areaPlaced / s.areaSheet) * 100.0;
        s.yMax = yMax(placed);
        s.elapsedMs = elapsedMs;
        return s;
    }


    private static double placedArea(List<PlacedPart> placed) {
        if (placed == null || placed.isEmpty()) return 0.0;
        List<Geometry> gs = new ArrayList<>(placed.size());
        for (PlacedPart p : placed) {
            if (p != null && p.polyPlaced != null && !p.polyPlaced.isEmpty()) {
                gs.add(p.polyPlaced);
            }
        }
        if (gs.isEmpty()) return 0.0;
        Geometry u = UnaryUnionOp.union(gs);
        return (u != null) ? u.getArea() : 0.0;
    }

    private static double yMax(List<PlacedPart> placed) {
        if (placed == null || placed.isEmpty()) return 0.0;
        double ymax = Double.NEGATIVE_INFINITY;
        for (PlacedPart p : placed) {
            if (p == null || p.polyPlaced == null || p.polyPlaced.isEmpty()) continue;
            Envelope env = p.polyPlaced.getEnvelopeInternal();
            if (env != null) ymax = Math.max(ymax, env.getMaxY());
        }
        return (ymax == Double.NEGATIVE_INFINITY) ? 0.0 : ymax;
    }

    /** -1 si a es mejor que b; +1 si peor; 0 empate. */
    private static int compare(Stats a, Stats b) {
        if (a == null && b == null) return 0;
        if (a == null) return 1;
        if (b == null) return -1;

        if (Double.compare(a.areaPlaced, b.areaPlaced) != 0) return (a.areaPlaced > b.areaPlaced) ? -1 : 1;
        if (a.placedCount != b.placedCount) return (a.placedCount > b.placedCount) ? -1 : 1;  
        if (Double.compare(a.yMax, b.yMax) != 0) return (a.yMax < b.yMax) ? -1 : 1;
        return 0;
    }

    private static Stats bestOf(List<Stats> list) {
        Stats best = null;
        for (Stats s : list) if (compare(s, best) < 0) best = s;
        return best;
    }

    // ===================== Saneamiento de soluciones (no solapes) =====================

    private static List<PlacedPart> sanitizeSolution(Polygon sheet, List<PlacedPart> sol) {
        if (sol == null || sol.isEmpty()) return Collections.emptyList();
        List<PlacedPart> out = new ArrayList<>();
        Geometry acc = null;
        for (PlacedPart p : sol) {
            if (p == null || p.polyPlaced == null || p.polyPlaced.isEmpty()) continue;
            if (!p.polyPlaced.within(sheet)) continue;
            if (acc == null) {
                acc = p.polyPlaced;
                out.add(clonePlaced(p));
            } else if (!p.polyPlaced.intersects(acc)) {
                out.add(clonePlaced(p));
                acc = UnaryUnionOp.union(Arrays.asList(acc, p.polyPlaced));
            }
        }
        return out;
    }

    private static PlacedPart clonePlaced(PlacedPart p) {
        PlacedPart q = new PlacedPart();
        q.id = p.id;
        q.polyPlaced = (p.polyPlaced != null) ? (Geometry) p.polyPlaced.copy() : null;
        q.angleDeg = p.angleDeg;
        q.mirrored = p.mirrored;
        q.dx = p.dx;
        q.dy = p.dy;
        return q;
    }

    // ===================== Logging helpers =====================

    private static PrintStream openLog(String logFile) throws Exception {
        if (logFile == null || logFile.isEmpty()) {
            new File(LOG_DIR).mkdirs();
            logFile = LOG_FILE;
        } else {
            File lf = new File(logFile).getAbsoluteFile();
            File parent = lf.getParentFile();
            if (parent != null) parent.mkdirs();
        }
        // PrintStream(OutputStream out, boolean autoFlush, String encoding)
        return new PrintStream(new FileOutputStream(logFile, true), true, StandardCharsets.UTF_8.name());
    }

    private static void head(PrintStream log, String title) {
        dash(log);
        line(log, title);
    }

    private static void summary(PrintStream log, String title) {
        dash(log);
        line(log, title);
    }

    private static void printStats(PrintStream log, Stats s) {
        dash(log);
        line(log, "RESUMEN (" + s.algo + ")");
        if (s.seed != 0) line(log, "Seed               : " + s.seed);
        line(log, "Piezas a posicionar: " + s.toPlace);
        line(log, "Piezas colocadas   : " + s.placedCount);
        line(log, "% aprovechamiento  : " + DF3.format(s.utilizationPct) + " %");
        line(log, "Área colocada      : " + DF3.format(s.areaPlaced));
        line(log, "Área plancha       : " + DF3.format(s.areaSheet));
        line(log, "Área desperdiciada : " + DF3.format(s.areaWaste));
        line(log, "Ymax               : " + DF3.format(s.yMax));
        line(log, "Tiempo " + s.algo + "     : " + s.elapsedMs + " ms (" + DF3.format(s.elapsedMs / 1000.0) + " s)");
        dash(log);
    }

    private static void renderSummary(PrintStream log, int partsCount, Stats best, long totalMs) {
        line(log, "Seed               : " + best.seed);
        line(log, "Piezas a posicionar: " + partsCount);
        line(log, "Piezas colocadas   : " + best.placedCount);
        line(log, "% aprovechamiento  : " + DF3.format(best.utilizationPct) + " %");
        line(log, "Área colocada      : " + DF3.format(best.areaPlaced));
        line(log, "Área plancha       : " + DF3.format(best.areaSheet));
        line(log, "Área desperdiciada : " + DF3.format(best.areaWaste));
        line(log, "Ymax               : " + DF3.format(best.yMax));
        line(log, "Tiempo total       : " + totalMs + " ms (" + DF3.format(totalMs / 1000.0) + " s)");
    }

    private static void dash(PrintStream log) {
        line(log, "------------------------------");
    }

    private static void line(PrintStream log, String msg) {
        log.println(TS.format(new Date()) + " | " + msg);
    }

    // ===================== POJO de métricas =====================
    private static class Stats {
        String algo;
        long   seed;
        List<PlacedPart> solution;
        int    toPlace;
        int    placedCount;
        double areaPlaced;
        double areaSheet;
        double areaWaste;
        double yMax;
        double utilizationPct;
        long   elapsedMs;
    }
}
