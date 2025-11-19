package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import java.io.*;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.*;

/**
 * Script de calibración automatizado para el algoritmo híbrido.
 * Permite ejecutar experimentos por partes para evitar tiempos de ejecución excesivos.
 */
public class CalibrationRunner {
    
    private static final DecimalFormat DF3 = new DecimalFormat("#,##0.000");
    private static final String RESULTS_DIR = "calibration_results";
    private static final String CONFIG_FILE = "calibration_config.csv";
    private static final String RESULTS_FILE = RESULTS_DIR + "/calibration_results.csv";
    
    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            System.out.println("Uso: CalibrationRunner <comando> [opciones]");
            System.out.println("Comandos:");
            System.out.println("  run <inicio> <fin>     - Ejecuta experimentos del rango [inicio, fin]");
            System.out.println("  run-all                 - Ejecuta todos los experimentos");
            System.out.println("  process                 - Procesa resultados y genera tablas");
            System.out.println("  list                    - Lista todos los experimentos configurados");
            return;
        }
        
        String command = args[0];
        new File(RESULTS_DIR).mkdirs();
        
        switch (command) {
            case "run":
                if (args.length < 3) {
                    System.out.println("Error: run requiere inicio y fin");
                    System.out.println("Ejemplo: CalibrationRunner run 0 9");
                    return;
                }
                int start = Integer.parseInt(args[1]);
                int end = Integer.parseInt(args[2]);
                runExperiments(start, end);
                break;
            case "run-all":
                runAllExperiments();
                break;
            case "process":
                processResults();
                break;
            case "list":
                listExperiments();
                break;
            default:
                System.out.println("Comando desconocido: " + command);
        }
    }
    
    private static void listExperiments() throws Exception {
        List<ExperimentConfig> configs = loadConfig();
        System.out.println("Total de experimentos configurados: " + configs.size());
        System.out.println("\nID | Descripción");
        System.out.println("---|------------");
        for (int i = 0; i < configs.size(); i++) {
            ExperimentConfig exp = configs.get(i);
            System.out.println(String.format("%3d | %s", i, exp.description));
        }
    }
    
    private static void runAllExperiments() throws Exception {
        List<ExperimentConfig> configs = loadConfig();
        runExperiments(0, configs.size() - 1);
    }
    
    private static void runExperiments(int start, int end) throws Exception {
        List<ExperimentConfig> configs = loadConfig();
        
        if (start < 0 || end >= configs.size() || start > end) {
            System.out.println("Error: Rango inválido. Debe estar entre 0 y " + (configs.size() - 1));
            return;
        }
        
        System.out.println("=== Ejecutando experimentos " + start + " a " + end + " ===");
        System.out.println("Total: " + (end - start + 1) + " experimentos");
        
        // Cargar datos de entrada
        String inputCsv = "input/dataset_L.csv";
        List<VertexRow> rows = CsvGeometryReader.readCsv(inputCsv);
        GeometryFactory gf = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING), 0);
        Polygon sheet = GeomUtils.buildSheet(rows, gf);
        if (sheet == null) throw new IllegalStateException("No se encontró SHEET");
        
        List<PartSpec> parts = GeomUtils.buildParts(rows, gf);
        if (parts.isEmpty()) throw new IllegalStateException("No se encontraron PART");
        
        double globalClearance = parts.stream().mapToDouble(p -> p.clearance).min().orElse(0.0);
        Polygon sheetInset = GeomUtils.insetSheet(sheet, globalClearance);
        
        List<ResultRow> results = new ArrayList<>();
        
        for (int expId = start; expId <= end; expId++) {
            ExperimentConfig exp = configs.get(expId);
            System.out.println("\n[" + (expId - start + 1) + "/" + (end - start + 1) + "] " + exp.description);
            
            for (int rep = 0; rep < exp.repetitions; rep++) {
                System.out.print("  Repetición " + (rep + 1) + "/" + exp.repetitions + "... ");
                
                long seed = System.nanoTime() + rep;
                long t0 = System.currentTimeMillis();
                
                String logFile = RESULTS_DIR + "/exp_" + expId + "_rep_" + rep + ".log";
                
                try {
                    List<PlacedPart> solution = HybridGraspMemeticPlacer.place(
                        sheetInset,
                        parts,
                        exp.graspSeeds,
                        exp.memeticTimeMs,
                        exp.graspAlpha,
                        logFile
                    );
                    
                    long t1 = System.currentTimeMillis();
                    long elapsedMs = t1 - t0;
                    
                    // Calcular métricas
                    double areaPlaced = calculateArea(solution);
                    double yMax = calculateYMax(solution);
                    int placedCount = solution != null ? solution.size() : 0;
                    double areaSheet = sheetInset.getArea();
                    double fitness = calculateFitness(areaPlaced, placedCount, yMax);
                    
                    ResultRow result = new ResultRow();
                    result.expId = expId;
                    result.rep = rep;
                    result.description = exp.description;
                    result.graspSeeds = exp.graspSeeds;
                    result.graspAlpha = exp.graspAlpha;
                    result.genMax = exp.genMax;
                    result.popSize = exp.popSize;
                    result.eliteSize = exp.eliteSize;
                    result.lsPasses = exp.lsPasses;
                    result.mutRotMax = exp.mutRotMax;
                    result.memeticTimeMs = exp.memeticTimeMs;
                    result.fitness = fitness;
                    result.areaPlaced = areaPlaced;
                    result.placedCount = placedCount;
                    result.yMax = yMax;
                    result.elapsedMs = elapsedMs;
                    
                    results.add(result);
                    System.out.println("OK (fitness=" + DF3.format(fitness) + ", tiempo=" + (elapsedMs/1000.0) + "s)");
                    
                } catch (Exception e) {
                    System.out.println("ERROR: " + e.getMessage());
                    e.printStackTrace();
                }
            }
        }
        
        // Guardar resultados
        saveResults(results, start, end);
        System.out.println("\n=== Resultados guardados en " + RESULTS_FILE + " ===");
    }
    
    private static void processResults() throws Exception {
        System.out.println("=== Procesando resultados ===");
        
        List<ResultRow> allResults = loadResults();
        if (allResults.isEmpty()) {
            System.out.println("No se encontraron resultados. Ejecuta primero algunos experimentos.");
            return;
        }
        
        // Agrupar por experimento
        Map<Integer, List<ResultRow>> byExp = new HashMap<>();
        for (ResultRow r : allResults) {
            byExp.computeIfAbsent(r.expId, k -> new ArrayList<>()).add(r);
        }
        
        // Calcular estadísticas por experimento
        List<ExperimentStats> stats = new ArrayList<>();
        for (Map.Entry<Integer, List<ResultRow>> entry : byExp.entrySet()) {
            List<ResultRow> expResults = entry.getValue();
            if (expResults.isEmpty()) continue;
            
            ExperimentStats stat = new ExperimentStats();
            stat.expId = entry.getKey();
            stat.description = expResults.get(0).description;
            stat.graspSeeds = expResults.get(0).graspSeeds;
            stat.graspAlpha = expResults.get(0).graspAlpha;
            stat.genMax = expResults.get(0).genMax;
            stat.popSize = expResults.get(0).popSize;
            stat.eliteSize = expResults.get(0).eliteSize;
            stat.lsPasses = expResults.get(0).lsPasses;
            stat.mutRotMax = expResults.get(0).mutRotMax;
            stat.memeticTimeMs = expResults.get(0).memeticTimeMs;
            
            // Calcular promedios
            stat.avgFitness = expResults.stream().mapToDouble(r -> r.fitness).average().orElse(0);
            stat.bestFitness = expResults.stream().mapToDouble(r -> r.fitness).max().orElse(0);
            stat.avgTime = expResults.stream().mapToLong(r -> r.elapsedMs).average().orElse(0) / 1000.0;
            stat.repetitions = expResults.size();
            
            stats.add(stat);
        }
        
        stats.sort(Comparator.comparingInt(s -> s.expId));
        
        // Generar tablas
        generateTables(stats);
        
        System.out.println("\n=== Tablas generadas en " + RESULTS_DIR + "/ ===");
    }
    
    private static void generateTables(List<ExperimentStats> stats) throws Exception {
        // Separar por tipo de experimento
        List<ExperimentStats> graspSeeds = new ArrayList<>();
        List<ExperimentStats> graspAlpha = new ArrayList<>();
        List<ExperimentStats> memetic = new ArrayList<>();
        List<ExperimentStats> memeticTime = new ArrayList<>();
        
        for (ExperimentStats s : stats) {
            if (s.description.contains("GRASP_SEEDS")) {
                graspSeeds.add(s);
            } else if (s.description.contains("GRASP_ALPHA")) {
                graspAlpha.add(s);
            } else if (s.description.contains("MEMETIC_TIME")) {
                memeticTime.add(s);
            } else if (s.description.contains("MEMETIC")) {
                memetic.add(s);
            }
        }
        
        // Tabla 1: Calibración GRASP_SEEDS
        if (!graspSeeds.isEmpty()) {
            writeTable("tabla_grasp_seeds.md", "Calibración GRASP_SEEDS", graspSeeds,
                new String[]{"GRASP_SEEDS", "α (fijo)", "Fitness Promedio", "Mejor Fitness", "Tiempo (s)"},
                s -> String.format("%d | %.1f | %.1f | %.1f | %.1f",
                    s.graspSeeds, s.graspAlpha, s.avgFitness, s.bestFitness, s.avgTime));
        }
        
        // Tabla 2: Calibración α
        if (!graspAlpha.isEmpty()) {
            writeTable("tabla_grasp_alpha.md", "Calibración α", graspAlpha,
                new String[]{"GRASP_SEEDS (fijo)", "α", "Fitness Promedio", "Mejor Fitness", "Tiempo (s)"},
                s -> String.format("%d | %.1f | %.1f | %.1f | %.1f",
                    s.graspSeeds, s.graspAlpha, s.avgFitness, s.bestFitness, s.avgTime));
        }
        
        // Tabla 3: Calibración Memético
        if (!memetic.isEmpty()) {
            writeTable("tabla_memetic.md", "Calibración Memético", memetic,
                new String[]{"GEN_MAX", "POP_SIZE", "ELITE_SIZE", "LS_PASSES", "MUT_ROT_MAX", "Fitness Promedio", "Mejor Fitness", "Tiempo (s)"},
                s -> String.format("%d | %d | %d | %d | %.1f | %.1f | %.1f | %.1f",
                    s.genMax, s.popSize, s.eliteSize, s.lsPasses, s.mutRotMax, s.avgFitness, s.bestFitness, s.avgTime));
        }
        
        // Tabla 4: Calibración Tiempo
        if (!memeticTime.isEmpty()) {
            writeTable("tabla_memetic_time.md", "Calibración MEMETIC_TIME_MS", memeticTime,
                new String[]{"Tiempo máximo (ms)", "Fitness Promedio", "Mejor Fitness", "Tiempo real (s)"},
                s -> String.format("%d | %.1f | %.1f | %.1f",
                    s.memeticTimeMs, s.avgFitness, s.bestFitness, s.avgTime));
        }
    }
    
    private static void writeTable(String filename, String title, List<ExperimentStats> stats,
                                   String[] headers, java.util.function.Function<ExperimentStats, String> formatter) throws Exception {
        File file = new File(RESULTS_DIR + "/" + filename);
        try (PrintWriter pw = new PrintWriter(new FileWriter(file, StandardCharsets.UTF_8))) {
            pw.println("# " + title);
            pw.println();
            pw.println("| " + String.join(" | ", headers) + " |");
            pw.println("|" + "---|".repeat(headers.length));
            for (ExperimentStats s : stats) {
                pw.println("| " + formatter.apply(s) + " |");
            }
        }
    }
    
    private static List<ExperimentConfig> loadConfig() throws Exception {
        List<ExperimentConfig> configs = new ArrayList<>();
        
        // Si no existe el archivo, crear configuración por defecto
        File configFile = new File(CONFIG_FILE);
        if (!configFile.exists()) {
            createDefaultConfig();
        }
        
        try (BufferedReader br = new BufferedReader(new FileReader(configFile, StandardCharsets.UTF_8))) {
            String line = br.readLine(); // Header
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split(",");
                if (parts.length < 11) continue;
                
                ExperimentConfig exp = new ExperimentConfig();
                exp.description = parts[0].trim();
                exp.repetitions = Integer.parseInt(parts[1].trim());
                exp.graspSeeds = Integer.parseInt(parts[2].trim());
                exp.graspAlpha = Double.parseDouble(parts[3].trim());
                exp.genMax = Integer.parseInt(parts[4].trim());
                exp.popSize = Integer.parseInt(parts[5].trim());
                exp.eliteSize = Integer.parseInt(parts[6].trim());
                exp.lsPasses = Integer.parseInt(parts[7].trim());
                exp.mutRotMax = Double.parseDouble(parts[8].trim());
                exp.memeticTimeMs = Long.parseLong(parts[9].trim());
                exp.iterMax = Integer.parseInt(parts[10].trim());
                exp.maxMillisTotal = Long.parseLong(parts[11].trim());
                
                configs.add(exp);
            }
        }
        
        return configs;
    }
    
    private static void createDefaultConfig() throws Exception {
        System.out.println("Creando archivo de configuración por defecto: " + CONFIG_FILE);
        
        try (PrintWriter pw = new PrintWriter(new FileWriter(CONFIG_FILE, StandardCharsets.UTF_8))) {
            pw.println("Descripción,Repeticiones,GRASP_SEEDS,GRASP_ALPHA,GEN_MAX,POP_SIZE,ELITE_SIZE,LS_PASSES,MUT_ROT_MAX,MEMETIC_TIME_MS,ITER_MAX,MAX_MILLIS_TOTAL");
            
            // Experimento 1: Variar GRASP_SEEDS (α fijo = 0.4)
            for (int seeds : new int[]{3, 5, 8, 10}) {
                pw.println(String.format("GRASP_SEEDS=%d (alpha=0.4),10,%d,0.4,100,36,4,2,2.0,120000,30,120000", seeds, seeds));
            }
            
            // Experimento 2: Variar α (GRASP_SEEDS fijo = 8)
            for (double alpha : new double[]{-1.0, 0.2, 0.4, 0.6, 0.8}) {
                pw.println(String.format("GRASP_ALPHA=%.1f (seeds=8),10,8,%.1f,100,36,4,2,2.0,120000,30,120000", alpha, alpha));
            }
            
            // Experimento 3: Variar GEN_MAX
            for (int genMax : new int[]{50, 100, 150, 200}) {
                pw.println(String.format("MEMETIC GEN_MAX=%d,10,8,0.4,%d,36,4,2,2.0,120000,30,120000", genMax, genMax));
            }
            
            // Experimento 4: Variar POP_SIZE
            for (int popSize : new int[]{24, 36, 50}) {
                pw.println(String.format("MEMETIC POP_SIZE=%d,10,8,0.4,100,%d,4,2,2.0,120000,30,120000", popSize, popSize));
            }
            
            // Experimento 5: Variar ELITE_SIZE
            for (int elite : new int[]{3, 4, 5}) {
                pw.println(String.format("MEMETIC ELITE_SIZE=%d,10,8,0.4,100,36,%d,2,2.0,120000,30,120000", elite, elite));
            }
            
            // Experimento 6: Variar LS_PASSES
            for (int lsPasses : new int[]{1, 2, 3}) {
                pw.println(String.format("MEMETIC LS_PASSES=%d,10,8,0.4,100,36,4,%d,2.0,120000,30,120000", lsPasses, lsPasses));
            }
            
            // Experimento 7: Variar MUT_ROT_MAX_DEG
            for (double rot : new double[]{2.0, 3.0}) {
                pw.println(String.format("MEMETIC MUT_ROT_MAX=%.1f,10,8,0.4,100,36,4,2,%.1f,120000,30,120000", rot, rot));
            }
            
            // Experimento 8: Variar MEMETIC_TIME_MS
            for (long time : new long[]{60000, 120000, 180000}) {
                pw.println(String.format("MEMETIC_TIME_MS=%d,10,8,0.4,100,36,4,2,2.0,%d,30,120000", time, time));
            }
        }
        
        System.out.println("Archivo de configuración creado. Edítalo según tus necesidades.");
    }
    
    private static void saveResults(List<ResultRow> results, int start, int end) throws Exception {
        File file = new File(RESULTS_FILE);
        boolean append = file.exists();
        
        try (PrintWriter pw = new PrintWriter(new FileWriter(file, StandardCharsets.UTF_8, append))) {
            if (!append) {
                pw.println("ExpID,Rep,Descripción,GRASP_SEEDS,GRASP_ALPHA,GEN_MAX,POP_SIZE,ELITE_SIZE,LS_PASSES,MUT_ROT_MAX,MEMETIC_TIME_MS,Fitness,AreaPlaced,PlacedCount,YMax,ElapsedMs");
            }
            
            for (ResultRow r : results) {
                pw.println(String.format("%d,%d,%s,%d,%.2f,%d,%d,%d,%d,%.1f,%d,%.2f,%.2f,%d,%.2f,%d",
                    r.expId, r.rep, r.description, r.graspSeeds, r.graspAlpha, r.genMax, r.popSize,
                    r.eliteSize, r.lsPasses, r.mutRotMax, r.memeticTimeMs, r.fitness, r.areaPlaced,
                    r.placedCount, r.yMax, r.elapsedMs));
            }
        }
    }
    
    private static List<ResultRow> loadResults() throws Exception {
        List<ResultRow> results = new ArrayList<>();
        File file = new File(RESULTS_FILE);
        if (!file.exists()) return results;
        
        try (BufferedReader br = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
            String line = br.readLine(); // Header
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] parts = line.split(",");
                if (parts.length < 16) continue;
                
                ResultRow r = new ResultRow();
                r.expId = Integer.parseInt(parts[0]);
                r.rep = Integer.parseInt(parts[1]);
                r.description = parts[2];
                r.graspSeeds = Integer.parseInt(parts[3]);
                r.graspAlpha = Double.parseDouble(parts[4]);
                r.genMax = Integer.parseInt(parts[5]);
                r.popSize = Integer.parseInt(parts[6]);
                r.eliteSize = Integer.parseInt(parts[7]);
                r.lsPasses = Integer.parseInt(parts[8]);
                r.mutRotMax = Double.parseDouble(parts[9]);
                r.memeticTimeMs = Long.parseLong(parts[10]);
                r.fitness = Double.parseDouble(parts[11]);
                r.areaPlaced = Double.parseDouble(parts[12]);
                r.placedCount = Integer.parseInt(parts[13]);
                r.yMax = Double.parseDouble(parts[14]);
                r.elapsedMs = Long.parseLong(parts[15]);
                
                results.add(r);
            }
        }
        
        return results;
    }
    
    // Funciones auxiliares para calcular métricas
    private static double calculateArea(List<PlacedPart> placed) {
        if (placed == null || placed.isEmpty()) return 0.0;
        double sum = 0.0;
        for (PlacedPart p : placed) {
            if (p != null && p.polyPlaced != null && !p.polyPlaced.isEmpty()) {
                sum += p.polyPlaced.getArea();
            }
        }
        return sum;
    }
    
    private static double calculateYMax(List<PlacedPart> placed) {
        if (placed == null || placed.isEmpty()) return 0.0;
        double ymax = Double.NEGATIVE_INFINITY;
        for (PlacedPart p : placed) {
            if (p != null && p.polyPlaced != null && !p.polyPlaced.isEmpty()) {
                Envelope env = p.polyPlaced.getEnvelopeInternal();
                if (env != null) ymax = Math.max(ymax, env.getMaxY());
            }
        }
        return (ymax == Double.NEGATIVE_INFINITY) ? 0.0 : ymax;
    }
    
    /**
     * Calcula un valor auxiliar para almacenamiento/visualización en CSV.
     * 
     * NOTA: La optimización real usa comparación lexicográfica pura (ver HybridMemetic.compare()),
     * NO este valor compuesto. Este método solo genera un número para facilitar el almacenamiento
     * y análisis estadístico en los archivos CSV, pero NO se usa para decidir qué solución es mejor.
     * 
     * La función objetivo es lexicográfica:
     * 1. Maximizar área colocada (prioridad 1)
     * 2. Maximizar cantidad de piezas colocadas (prioridad 2, solo si área es igual)
     * 3. Minimizar YMax (prioridad 3, solo si área y cantidad son iguales)
     */
    private static double calculateFitness(double areaPlaced, int placedCount, double yMax) {
        // Valor auxiliar aproximado para CSV (NO usado en optimización)
        return areaPlaced * 1000 + placedCount * 10 - yMax * 0.01;
    }
    
    // Clases auxiliares
    static class ExperimentConfig {
        String description;
        int repetitions;
        int graspSeeds;
        double graspAlpha;
        int genMax;
        int popSize;
        int eliteSize;
        int lsPasses;
        double mutRotMax;
        long memeticTimeMs;
        int iterMax;
        long maxMillisTotal;
    }
    
    static class ResultRow {
        int expId;
        int rep;
        String description;
        int graspSeeds;
        double graspAlpha;
        int genMax;
        int popSize;
        int eliteSize;
        int lsPasses;
        double mutRotMax;
        long memeticTimeMs;
        double fitness;
        double areaPlaced;
        int placedCount;
        double yMax;
        long elapsedMs;
    }
    
    static class ExperimentStats {
        int expId;
        String description;
        int graspSeeds;
        double graspAlpha;
        int genMax;
        int popSize;
        int eliteSize;
        int lsPasses;
        double mutRotMax;
        long memeticTimeMs;
        double avgFitness;
        double bestFitness;
        double avgTime;
        int repetitions;
    }
}




