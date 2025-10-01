package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.operation.union.UnaryUnionOp;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;
//PARA ESTE ALGORITMO SE USAN HYBRIDGRASP y HYBRIDMEMETIC
/**
 * ===================== HybridGraspMemeticPlacer =====================
 * Este algoritmo híbrido combina GRASP y Memético en varias fases:
 *
 * 1) GRASP Multistart: se ejecuta varias veces GRASP con diferentes seeds,
 *    generando múltiples layouts iniciales factibles.
 * 2) Selección de la mejor seed: se toma el layout que más piezas colocó.
 * 3) Algoritmo Memético: esa mejor seed se pasa como población inicial
 *    al memético (similar a un genético con búsqueda local).
 * 4) Selección base: se queda con la mejor solución entre GRASP y Memético.
 * 5) Fusión (opcional): intenta insertar piezas adicionales en huecos libres
 *    usando lógica de GRASP sobre la solución base. (que sea opcional se 
 *    refiere a que siempre se hace pero no siempre la solución mejora).
 *
 * Resultado: GRASP aporta diversidad y factibilidad, el Memético explora y mejora,
 * y la fusión agrega un último refinamiento para aumentar piezas colocadas.
 */
public class HybridGraspMemeticPlacer {

    // ====== Parámetros del híbrido ======
    private static final int    GRASP_SEEDS      = 5;     // número de seeds GRASP
    private static final long   BASE_SEED        = System.nanoTime();
    private static final double FUSION_MAX_FRAC  = 0.12;  // fracción máxima de piezas candidatas en fusión
    private static final int    FUSION_MAX_TRIES = 1;     // intentos de inserción de fusión

    private static final String LOG_FILE         = "LogHybrid.txt";
    private static final SimpleDateFormat TS     = new SimpleDateFormat("HH:mm:ss.SSS");
    private static final DecimalFormat DF        = new DecimalFormat("#,##0.###");

    private static PrintStream tee;

    /** === Punto de entrada principal: orquesta GRASP + Memético + Fusión === */
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs) {
        initLogger(false);
        long t0 = System.currentTimeMillis();
        log("==== HYB GRASP+MEMETIC | START ====");
        log("CFG | GRASP_SEEDS=%d  FUSION_MAX_FRAC=%.2f  FUSION_MAX_TRIES=%d",
                GRASP_SEEDS, FUSION_MAX_FRAC, FUSION_MAX_TRIES);

        List<PlacedPart> finalLayout = Collections.emptyList();
        long chosenSeed = BASE_SEED; // para reporte
        String chosenAlgo = "HYBRID"; // etiqueta
        try {
            // ---------- 1) MULTISTART GRASP ----------
            // Genera varias soluciones con GRASP usando seeds distintos
            long tGrasp = System.currentTimeMillis();
            List<List<PlacedPart>> seeds = HybridGrasp.multiStart(sheetInset, specs, GRASP_SEEDS, BASE_SEED, tee);

            // Se elige la mejor seed (la que más piezas colocó)
            int bestIdx = -1, bestPlaced = -1;
            for (int i = 0; i < seeds.size(); i++) {
                int placed = (seeds.get(i) == null ? 0 : seeds.get(i).size());
                if (placed > bestPlaced) { bestPlaced = placed; bestIdx = i; }
            }
            List<PlacedPart> bestSeed = (bestIdx >= 0 ? seeds.get(bestIdx) : Collections.emptyList());
            chosenSeed = BASE_SEED + Math.max(0, bestIdx);
            int placedSeed = (bestSeed == null ? 0 : bestSeed.size());
            log("GRASP | best-seed idx=%d seed=%d placed=%d  (%.3fs)",
                    bestIdx, chosenSeed, placedSeed, secsSince(tGrasp));

            // Si GRASP ya colocó todas las piezas, no hace falta el memético
            int totalParts = totalQuantity(specs);
            if (placedSeed >= totalParts && totalParts > 0) {
                log("SKIP MEMETIC | GRASP colocó todas las piezas (%d/%d).", placedSeed, totalParts);
                finalLayout = safe(bestSeed);
                chosenAlgo = "GRASP";
                printFinalSummary("HYBRID", sheetInset, specs, finalLayout, chosenSeed, t0);
                return finalLayout;
            }

            // ---------- 2) MEMÉTICO con SEMILLA ----------
            // Aquí se ejecuta el algoritmo memético, usando como semilla el mejor layout de GRASP
            long tMem = System.currentTimeMillis();
            HybridMemetic.Config mc = new HybridMemetic.Config();
            mc.pop = 12; mc.gen = 12;
            mc.budgetMs = 0L;
            mc.captureStdout = true;         // Redirige la salida estándar al log
            mc.postRepairTries = 0;
            mc.rngSeed = BASE_SEED + 999;

            log("MEM  | run with seed (si MemeticPlacer soporta seed, la usará)");
            List<PlacedPart> mem = HybridMemetic.run(sheetInset, specs, bestSeed, mc, tee, HybridGrasp::reinsertOnFree);
            int placedMem = (mem == null ? 0 : mem.size());
            log("MEM  | end placed=%d  (%.3fs)", placedMem, secsSince(tMem));

            // ---------- 3) SELECCIÓN BASE ----------
            // Se queda con la mejor solución entre GRASP y Memético
            List<PlacedPart> base = (placedMem > placedSeed ? mem : bestSeed);
            String chosenBase = (base == mem ? "MEMETIC" : "GRASP");
            int basePlaced = (base == null ? 0 : base.size());
            log("BASE | chosen=%s placed=%d", chosenBase, basePlaced);

            // ---------- 4) FUSIÓN OPCIONAL ----------
            // Intenta mejorar aún más insertando piezas con lógica GRASP en huecos libres
            List<PlacedPart> other = (base == mem ? bestSeed : mem);
            List<PlacedPart> fused = tryFusionInsert(sheetInset, specs, base, other);
            int fusedPlaced = (fused == null ? 0 : fused.size());
            if (fusedPlaced > basePlaced) {
                log("FUSE | improved: %d -> %d", basePlaced, fusedPlaced);
                finalLayout = fused;
                chosenAlgo  = "HYBRID"; // indica que se aplicó la fusión
            } else {
                log("FUSE | no improvement (kept %s)", chosenBase);
                finalLayout = base;
                chosenAlgo  = chosenBase;
            }

            // ---------- 5) RESUMEN FINAL ----------
            printFinalSummary("HYBRID", sheetInset, specs, finalLayout, chosenSeed, t0);
            return finalLayout;

        } catch (Throwable ex) {
            log("!! HYB ERROR | %s", ex.toString());
            log(stackTraceOf(ex));
            log("Fallback => retorno vacío.");
            printFinalSummary("HYBRID", sheetInset, specs, finalLayout, chosenSeed, t0);
            return Collections.emptyList();
        } finally {
            closeLogger();
        }
    }

    private static List<PlacedPart> tryFusionInsert(Polygon sheetInset,
                                                    List<PartSpec> allSpecs,
                                                    List<PlacedPart> base,
                                                    List<PlacedPart> other) {
        if (base == null) base = Collections.emptyList();
        if (other == null || other.isEmpty()) {
            log("FUSE | skip (other vacío)");
            return base;
        }

        // Construye conjunto de IDs de las piezas que ya están en base
        Set<Object> idsBase = new HashSet<>();
        for (PlacedPart pp : base) {
            Object id = partIdOf(pp);
            if (id != null) idsBase.add(id);
        }

        // Candidatas = piezas que están en "other" pero no en "base"
        List<PartSpec> candidates = new ArrayList<>();
        for (PlacedPart pp : other) {
            Object id = partIdOf(pp);
            if (id != null && !idsBase.contains(id)) {
                PartSpec spec = specById(allSpecs, id);
                if (spec != null) candidates.add(spec);
            }
        }
        if (candidates.isEmpty()) {
            log("FUSE | no candidates (ambas colocan lo mismo)");
            return base;
        }

        // Limita cuántas piezas candidatas puede intentar insertar
        int maxAdds = Math.max(1, (int) Math.round(FUSION_MAX_FRAC * candidates.size()));
        if (candidates.size() > maxAdds) candidates = candidates.subList(0, maxAdds);

        // Calcula espacio libre actual de la plancha
        Geometry free = safeFreeSpaceOfLayout(sheetInset, base);
        if (!(free instanceof Polygon) && !(free instanceof MultiPolygon)) {
            log("FUSE | no free-space polygonal (skip)");
            return base;
        }

        // Copia de la solución base, sobre la que se van insertando piezas extra
        List<PlacedPart> current = new ArrayList<>(base);
        long seed0 = BASE_SEED + 202;
        int inserted = 0;

        // Intenta insertar cada candidata de a una
        for (int i = 0; i < candidates.size(); i++) {
            PartSpec spec = candidates.get(i);
            Polygon freePoly = toPolygonOrHull(free);
            if (freePoly == null) {
                log("FUSE | free no-poligonal (stop)");
                break;
            }

            // Reintenta colocar una pieza con GRASP en el espacio libre
            List<PartSpec> singleton = Collections.singletonList(spec);
            List<PlacedPart> placedOne = HybridGrasp.reinsertOnFree(freePoly, singleton, seed0 + i, tee);

            // Si se pudo colocar, se agrega a la solución y se ajusta el espacio libre
            if (placedOne != null && !placedOne.isEmpty()) {
                PlacedPart q = placedOne.get(0);
                current.add(q);
                Geometry fp = footprintOf(q);
                if (fp != null) {
                    free = safeBuffer(safeDifference(free, safeBuffer(fp, +eps(freePoly))), 0);
                }
                inserted++;
                if (inserted >= FUSION_MAX_TRIES * candidates.size()) break;
            }
        }

        log("FUSE | inserted=%d", inserted);
        return current;
    }


    // Calcula el espacio libre de la plancha dado un layout actual
    private static Geometry safeFreeSpaceOfLayout(Polygon sheetInset, List<PlacedPart> layout) {
        try {
            List<Geometry> geoms = new ArrayList<>();
            for (PlacedPart pp : layout) {
                Geometry fp = footprintOf(pp);
                if (fp != null) geoms.add(fp);
            }
            Geometry occ = geoms.isEmpty() ? null : UnaryUnionOp.union(geoms);
            if (occ == null) return sheetInset;
            Geometry dilated = safeBuffer(occ, +eps(sheetInset)); // expande un poquito para evitar errores de borde
            return safeDifference(sheetInset, dilated);
        } catch (Throwable t) {
            log("FUSE | free-space error: %s", t.toString());
            return sheetInset;
        }
    }

    // Obtiene la geometría (footprint) de una pieza colocada
    private static Geometry footprintOf(PlacedPart pp) {
        if (pp == null) return null;
        try {
            return (Geometry) GeomUtils.class
                    .getMethod("footprint", pp.getClass())
                    .invoke(null, pp);
        } catch (NoSuchMethodException ignore) {
            try { return (Geometry) pp.getClass().getMethod("getGeometry").invoke(pp); }
            catch (Exception e2) {
                try { return (Geometry) pp.getClass().getMethod("getPolygon").invoke(pp); }
                catch (Exception e3) { return null; }
            }
        } catch (Exception e) {
            return null;
        }
    }

    // Obtiene el ID de una pieza colocada
    private static Object partIdOf(PlacedPart pp) {
        if (pp == null) return null;
        try { return pp.getClass().getMethod("getId").invoke(pp); }
        catch (Exception ignore) {
            try { return pp.getClass().getMethod("id").invoke(pp); }
            catch (Exception ignore2) { return null; }
        }
    }

    // Busca en la lista de specs la que corresponde al ID dado
    private static PartSpec specById(List<PartSpec> specs, Object id) {
        if (id == null || specs == null) return null;
        for (PartSpec ps : specs) {
            try {
                Object sid = ps.getClass().getMethod("getId").invoke(ps);
                if (id.equals(sid)) return ps;
            } catch (Exception ignore) {
                try {
                    Object sid2 = ps.getClass().getMethod("id").invoke(ps);
                    if (id.equals(sid2)) return ps;
                } catch (Exception ignore2) {}
            }
        }
        return null;
    }

    // Convierte una geometría en polígono (o su convex hull si no lo es)
    private static Polygon toPolygonOrHull(Geometry g) {
        if (g == null) return null;
        if (g instanceof Polygon) return (Polygon) g;
        Geometry hull = g.convexHull();
        return (hull instanceof Polygon) ? (Polygon) hull : null;
    }

    // Operación segura: diferencia entre geometrías
    private static Geometry safeDifference(Geometry a, Geometry b) {
        try { return a.difference(b); } catch (Throwable t) { return a; }
    }

    // Operación segura: buffer de geometría
    private static Geometry safeBuffer(Geometry a, double d) {
        try { return a.buffer(d); } catch (Throwable t) { return a; }
    }

    // Calcula tolerancia (epsilon) relativa al tamaño de la plancha
    private static double eps(Polygon sheet) {
        Envelope e = sheet.getEnvelopeInternal();
        double minDim = Math.min(e.getWidth(), e.getHeight());
        return Math.max(1e-7, minDim / 1e7);
    }

    private static int totalQuantity(List<PartSpec> specs) {
        if (specs == null) return 0;
        int sum = 0;
        for (PartSpec ps : specs) {
            try {
                // Intenta acceder vía getter
                Object q = ps.getClass().getMethod("getQty").invoke(ps);
                if (q instanceof Number) { sum += ((Number) q).intValue(); continue; }
            } catch (Exception ignore) {}
            try {
                // Intenta acceder como campo público
                Object qf = ps.getClass().getField("qty").get(ps);
                if (qf instanceof Number) { sum += ((Number) qf).intValue(); continue; }
            } catch (Exception ignore) {}
            // Si no se encuentra, asume 1 para no subreportar
            sum += 1;
        }
        return sum;
    }

    // Imprime el resumen final del híbrido con formato consistente
    private static void printFinalSummary(String label,
                                          Polygon sheetInset,
                                          List<PartSpec> specs,
                                          List<PlacedPart> layout,
                                          long seedValue,
                                          long tStartMs) {
        // Calcula métricas
        int piezasObjetivo = totalQuantity(specs);
        int piezasColocadas = (layout == null ? 0 : layout.size());

        double areaSheet = (sheetInset == null ? 0.0 : sheetInset.getArea());

        Geometry occ = unionFootprints(layout);
        double areaColocada = (occ == null ? 0.0 : safeArea(occ));
        double areaDesperdiciada = Math.max(0.0, areaSheet - areaColocada);
        double aprovechamiento = (areaSheet <= 0.0 ? 0.0 : (areaColocada / areaSheet) * 100.0);

        double yMax = 0.0;
        if (occ != null) {
            Envelope env = occ.getEnvelopeInternal();
            yMax = env.getMaxY();
        }

        long elapsedMs = System.currentTimeMillis() - tStartMs;
        double elapsedSec = elapsedMs / 1000.0;

        // Imprime estilo tabla con timestamp
        String ts = TS.format(new Date());
        log("------------------------------");
        log("%s RESUMEN FINAL (%s)", "[" + ts + "]", label.toUpperCase());
        log("%s Seed               : %s", "[" + ts + "]", String.valueOf(seedValue));
        log("%s Piezas a posicionar: %s", "[" + ts + "]", DF.format(piezasObjetivo));
        log("%s Piezas colocadas   : %s", "[" + ts + "]", DF.format(piezasColocadas));
        log("%s %% aprovechamiento  : %s %%", "[" + ts + "]", DF.format(aprovechamiento));
        log("%s Área colocada      : %s", "[" + ts + "]", DF.format(areaColocada));
        log("%s Área plancha       : %s", "[" + ts + "]", DF.format(areaSheet));
        log("%s Área desperdiciada : %s", "[" + ts + "]", DF.format(areaDesperdiciada));
        log("%s Ymax               : %s", "[" + ts + "]", DF.format(yMax));
        log("%s Tiempo total       : %s ms (%s s)", "[" + ts + "]", DF.format(elapsedMs), DF.format(elapsedSec));
        log("------------------------------");
    }

    // Une footprints de todas las piezas para calcular métricas de área
    private static Geometry unionFootprints(List<PlacedPart> layout) {
        if (layout == null || layout.isEmpty()) return null;
        List<Geometry> geoms = new ArrayList<>();
        for (PlacedPart pp : layout) {
            Geometry g = footprintOf(pp);
            if (g != null) geoms.add(g);
        }
        if (geoms.isEmpty()) return null;
        try { return UnaryUnionOp.union(geoms); }
        catch (Throwable t) { return null; }
    }

    // Obtiene área de forma segura
    private static double safeArea(Geometry g) {
        try { return g.getArea(); } catch (Throwable t) { return 0.0; }
    }

    /* ===================== Logging / util ===================== */

    // Devuelve lista segura (no nula)
    private static List<PlacedPart> safe(List<PlacedPart> v) {
        return (v == null ? Collections.emptyList() : v);
    }

    // Inicializa logger: duplica salida en consola y archivo
    private static void initLogger(boolean append) {
        try {
            PrintStream file = new PrintStream(new FileOutputStream(LOG_FILE, append), true, StandardCharsets.UTF_8);
            tee = new TeePrintStream(System.out, file);
        } catch (Exception e) {
            tee = System.out; // fallback: solo consola
        }
    }

    // Cierra logger
    private static void closeLogger() {
        try { if (tee != null) tee.flush(); } catch (Exception ignore) {}
        try { if (tee != null) tee.close(); } catch (Exception ignore) {}
    }

    // Imprime línea de log con timestamp
    private static void log(String fmt, Object... args) {
        String ts = TS.format(new Date());
        tee.println("[" + ts + "] " + String.format(fmt, args));
    }

    // Calcula segundos desde t0
    private static double secsSince(long t0) {
        return (System.currentTimeMillis() - t0) / 1000.0;
    }

    // Devuelve stacktrace como String
    private static String stackTraceOf(Throwable t) {
        try (java.io.StringWriter sw = new java.io.StringWriter();
             java.io.PrintWriter pw = new java.io.PrintWriter(sw)) {
            t.printStackTrace(pw);
            return sw.toString();
        } catch (Exception e) {
            return t.toString();
        }
    }

    static class TeePrintStream extends PrintStream {
        private final PrintStream a, b;

        TeePrintStream(PrintStream out, PrintStream file) {
            super(out);
            this.a = out; this.b = file;
        }

        @Override public void println(String x) {
            a.println(x); b.println(x);
        }

        @Override public void print(String s) {
            a.print(s); b.print(s);
        }

        @Override public PrintStream printf(String format, Object... args) {
            a.printf(format, args);
            b.printf(format, args);
            return this;
        }

        @Override public void flush() {
            try { a.flush(); } catch (Exception ignore) {}
            try { b.flush(); } catch (Exception ignore) {}
        }

        @Override public void close() {
            try { a.close(); } catch (Exception ignore) {}
            try { b.close(); } catch (Exception ignore) {}
        }
    }
}
