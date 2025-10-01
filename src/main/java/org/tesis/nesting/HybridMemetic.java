package org.tesis.nesting;

import org.locationtech.jts.geom.Polygon;

import java.io.PrintStream;
import java.util.*;


public final class HybridMemetic {

    private HybridMemetic() {}

    public static class Config {
        public int    pop            = 12;                 // tamaño de población
        public int    gen            = 12;                 // número de generaciones
        public long   budgetMs       = 0L;                 // 0 => sin tope desde este lado
        public long   rngSeed        = System.nanoTime();  // semilla base
        public boolean captureStdout = true;               // si true, redirijo System.out al log
        public int    postRepairTries = 0;                 // intentos de post-repair (0 => desactivado)
        public double postRepairFrac  = 0.10;              // % de specs candidato a reinsertar en post
    }


    @FunctionalInterface
    public interface RepairOperator {
        List<PlacedPart> repair(Polygon freeSpace, List<PartSpec> subset, long seed, PrintStream log);
    }

    public static List<PlacedPart> run(Polygon sheetInset,
                                       List<PartSpec> specs,
                                       List<PlacedPart> seedLayout,
                                       Config cfg,
                                       PrintStream log,
                                       RepairOperator repairOp) {
        Objects.requireNonNull(sheetInset, "sheetInset");
        Objects.requireNonNull(specs, "specs");

        // Si no me pasan config, uso una por defecto.
        if (cfg == null) cfg = new Config();

        // Si no me pasan operador de reparación, por defecto uso el de GRASP.
        if (repairOp == null) {
            repairOp = HybridGrasp::reinsertOnFree;
        }

        // Guardo stdout actual por si tengo que redirigir.
        PrintStream oldOut = System.out;
        try {
            // Si me piden capturar stdout, mando System.out al stream de log.
            if (cfg.captureStdout && log != null) System.setOut(log);

            // 1) Primero pruebo sobrecargas del memético que aceptan semilla.
            List<PlacedPart> ans = tryMemeticWithSeed(sheetInset, specs, seedLayout, cfg, log);
            if (ans != null) {
                // Si el memético devolvió algo y está activado el post-repair, hago un intento rápido.
                if (cfg.postRepairTries > 0) {
                    ans = postRepairCheap(sheetInset, specs, ans, cfg, repairOp, log);
                }
                return ans;
            }

            // 2) Si no encontré firma con semilla, intento la firma clásica.
            logMsg(log, "HYB_MEMETIC | fallback -> MemeticPlacer.place(sheet, specs)");
            ans = callMemeticClassic(sheetInset, specs, cfg, log);

            // Aplico el mismo post-repair opcional si corresponde.
            if (cfg.postRepairTries > 0) {
                ans = postRepairCheap(sheetInset, specs, ans, cfg, repairOp, log);
            }
            return (ans == null ? Collections.emptyList() : ans);

        } finally {
            // Restauro stdout si lo toqué.
            try { if (cfg.captureStdout) System.setOut(oldOut); } catch (Exception ignore) {}
        }
    }

    @SuppressWarnings("unchecked")
    private static List<PlacedPart> tryMemeticWithSeed(Polygon sheet,
                                                       List<PartSpec> specs,
                                                       List<PlacedPart> seedLayout,
                                                       Config cfg,
                                                       PrintStream log) {
        // Si no hay semilla no tiene sentido entrar acá.
        if (seedLayout == null || seedLayout.isEmpty()) {
            logMsg(log, "HYB_MEMETIC | seedLayout vacío, salto sobrecargas con seed.");
            return null;
        }

        try {
            // Firma completa: place(sheet, specs, cfg, seedLayout, budgetMs, log)
            return (List<PlacedPart>) MemeticPlacer.class
                    .getMethod("place", Polygon.class, List.class, Config.class, List.class, long.class, PrintStream.class)
                    .invoke(null, sheet, specs, cfg, seedLayout, cfg.budgetMs, log);
        } catch (NoSuchMethodException ignore) {
            // No pasa nada, pruebo la siguiente.
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
            return null;
        }

        try {
            // Variante sin log: place(sheet, specs, cfg, seedLayout, budgetMs)
            return (List<PlacedPart>) MemeticPlacer.class
                    .getMethod("place", Polygon.class, List.class, Config.class, List.class, long.class)
                    .invoke(null, sheet, specs, cfg, seedLayout, cfg.budgetMs);
        } catch (NoSuchMethodException ignore) {
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
            return null;
        }

        try {
            // Versión mínima con semilla: place(sheet, specs, seedLayout)
            return (List<PlacedPart>) MemeticPlacer.class
                    .getMethod("place", Polygon.class, List.class, List.class)
                    .invoke(null, sheet, specs, seedLayout);
        } catch (NoSuchMethodException ignore) {
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
            return null;
        }

        // Si llegué hasta acá, no encontré ninguna firma con seed.
        return null;
    }

    @SuppressWarnings("unchecked")
    private static List<PlacedPart> callMemeticClassic(Polygon sheet,
                                                       List<PartSpec> specs,
                                                       Config cfg,
                                                       PrintStream log) {
        try {
            // Intento una firma más completa: place(sheet, specs, cfg, budgetMs, log)
            return (List<PlacedPart>) MemeticPlacer.class
                    .getMethod("place", Polygon.class, List.class, Config.class, long.class, PrintStream.class)
                    .invoke(null, sheet, specs, cfg, cfg.budgetMs, log);
        } catch (NoSuchMethodException ignore) {
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
        }

        try {
            // Firma base: place(sheet, specs)
            return (List<PlacedPart>) MemeticPlacer.class
                    .getMethod("place", Polygon.class, List.class)
                    .invoke(null, sheet, specs);
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
            return Collections.emptyList();
        }
    }

    private static List<PlacedPart> postRepairCheap(Polygon sheet,
                                                    List<PartSpec> specs,
                                                    List<PlacedPart> layout,
                                                    Config cfg,
                                                    RepairOperator repairOp,
                                                    PrintStream log) {
        if (cfg.postRepairTries <= 0 || specs == null || specs.isEmpty()) return layout;

        logMsg(log, "HYB_MEMETIC | post-repair: tries=%d frac=%.2f", cfg.postRepairTries, cfg.postRepairFrac);

        // Selecciono un subconjunto pseudo-aleatorio de specs (tamaño ≈ postRepairFrac * |specs|).
        int want = Math.max(1, (int) Math.round(cfg.postRepairFrac * specs.size()));
        List<PartSpec> shuffled = new ArrayList<>(specs);
        Collections.shuffle(shuffled, new Random(cfg.rngSeed));
        List<PartSpec> subset = shuffled.subList(0, Math.min(want, shuffled.size()));

        // Arranco del layout actual; si está en null, lo trato como vacío.
        List<PlacedPart> best = (layout == null ? Collections.emptyList() : layout);
        int bestPlaced = best.size();
        long seed = cfg.rngSeed;

        // Hago unos pocos intentos. Si el repair devuelve algo, concateno y comparo.
        for (int t = 0; t < cfg.postRepairTries; t++) {
            List<PlacedPart> repaired = repairOp.repair(sheet, subset, seed + t, log);
            if (repaired != null && !repaired.isEmpty()) {
                List<PlacedPart> candidate = new ArrayList<>(best);
                candidate.addAll(repaired);
                if (candidate.size() > bestPlaced) {
                    best = candidate;
                    bestPlaced = candidate.size();
                    logMsg(log, "HYB_MEMETIC | post-repair sube a placed=%d", bestPlaced);
                }
            }
        }
        return best;
    }

    private static void logMsg(PrintStream log, String fmt, Object... args) {
        if (log != null) log.println(String.format(fmt, args));
    }
}
