package org.tesis.nesting;

import org.locationtech.jts.geom.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.util.List;


public class TesisNesting {
    
    // ===================== Configuración =====================
    // Número de veces que se ejecutará cada algoritmo
    private static final int NUM_ITERACIONES = 1;
    
    // Archivos CSV donde se guardarán los resultados
    private static final String RESULTADOS_HIBRIDOS_CSV = "out/resultados_hibridos.csv";
    private static final String RESULTADOS_GRASP_CSV = "out/resultados_grasp.csv";
    
    private static final DecimalFormat DF3 = new DecimalFormat("#,##0.000");
 
    public static void main(String[] args) throws Exception {
        //Acá se puede cambiar el nombre de el svg resultante.
        //String inputCsv  = "input/datosVariados.csv";
        String inputCsv  = "input/dataset_L.csv";
        String outputDir = "out";
        String outputSvg = outputDir + "/salidaHybridPlacer.svg";

        if (args.length >= 1) inputCsv = args[0];
        if (args.length >= 2) outputSvg = args[1];

        new File(outputDir).mkdirs();

        // 1) Leer CSV
        List<VertexRow> rows = CsvGeometryReader.readCsv(inputCsv);

        // 2) Reconstruir geometrías
        GeometryFactory gf = new GeometryFactory(new PrecisionModel(PrecisionModel.FLOATING), 0);
        Polygon sheet = GeomUtils.buildSheet(rows, gf);
        if (sheet == null) throw new IllegalStateException("No se encontró SHEET en " + inputCsv);

        List<PartSpec> parts = GeomUtils.buildParts(rows, gf);
        if (parts.isEmpty()) throw new IllegalStateException("No se encontraron PART en " + inputCsv);

        // 3) Inset de la plancha según clearance global mínimo
        double globalClearance = parts.stream().mapToDouble(p -> p.clearance).min().orElse(0.0);
        Polygon sheetInset = GeomUtils.insetSheet(sheet, globalClearance);

        // 4) Ejecutar algoritmo GRASP múltiples veces (solo GRASP, no híbrido)
        ejecutarAlgoritmoHibrido(sheetInset, sheet, parts, outputSvg);
    }
    
    private static void ejecutarAlgoritmoHibrido(Polygon sheetInset, Polygon sheet, List<PartSpec> parts, String outputSvg) throws Exception {
        // Inicializar archivo CSV (escribir encabezados si no existe)
        boolean archivoExiste = new File(RESULTADOS_HIBRIDOS_CSV).exists();
        try (PrintWriter csvWriter = new PrintWriter(
                new OutputStreamWriter(new FileOutputStream(RESULTADOS_HIBRIDOS_CSV, true), StandardCharsets.UTF_8))) {
            
            // Escribir encabezados solo si el archivo no existe
            if (!archivoExiste) {
                csvWriter.println("Corrida,Seed,Algoritmo,Piezas_a_posicionar,Piezas_colocadas,Area_colocada,Area_plancha,Area_desperdiciada,Ymax,Utilizacion_%,Tiempo_ms,Tiempo_s");
            }
            
            // Ejecutar algoritmo híbrido múltiples veces
            System.out.println("========================================");
            System.out.println("Iniciando " + NUM_ITERACIONES + " ejecuciones del algoritmo HÍBRIDO");
            System.out.println("========================================");
            
            HybridResult ultimoResult = null;
            
            for (int corrida = 1; corrida <= NUM_ITERACIONES; corrida++) {
                System.out.println("\n>>> HÍBRIDO - CORRIDA " + corrida + " de " + NUM_ITERACIONES + " <<<");
                
                // Ejecutar algoritmo híbrido y obtener estadísticas
                HybridResult result = HybridGraspMemeticPlacer.placeWithStats(sheetInset, parts);
                result.corrida = corrida;
                ultimoResult = result; // Guardar para el SVG final
                
                // Guardar resultado en CSV
                escribirResultadoCSV(csvWriter, result);
                csvWriter.flush(); // Asegurar que se escriba inmediatamente
                
                System.out.println(">>> Resultado guardado en CSV: Corrida " + corrida + 
                                 " | Colocadas: " + result.placedCount + 
                                 " | Area: " + DF3.format(result.areaPlaced) + 
                                 " | Ymax: " + DF3.format(result.yMax) + 
                                 " | Utilización: " + DF3.format(result.utilizationPct) + "%");
            }
            
            System.out.println("\n========================================");
            System.out.println("Todas las corridas HÍBRIDAS completadas");
            System.out.println("Resultados guardados en: " + RESULTADOS_HIBRIDOS_CSV);
            System.out.println("========================================");
            
            // Exportar SVG de la última corrida (opcional, para visualización)
            if (ultimoResult != null) {
                String svg = SvgWriter.toSVG(sheet, ultimoResult.solution);
                try (Writer w = new OutputStreamWriter(new FileOutputStream(outputSvg), StandardCharsets.UTF_8)) {
                    w.write(svg);
                }
                System.err.println("SVG generado en: " + outputSvg);
            }
        }
    }
    
    private static void ejecutarAlgoritmoGrasp(Polygon sheetInset, Polygon sheet, List<PartSpec> parts, String outputSvg) throws Exception {
        // Verificar si el archivo existe y contar corridas existentes
        File archivoCSV = new File(RESULTADOS_GRASP_CSV);
        boolean archivoExiste = archivoCSV.exists();
        int corridaInicial = 1;
        
        if (archivoExiste) {
            // Contar cuántas corridas ya existen en el archivo
            try (BufferedReader br = new BufferedReader(
                    new InputStreamReader(new FileInputStream(archivoCSV), StandardCharsets.UTF_8))) {
                String header = br.readLine(); // Leer encabezado
                if (header != null) {
                    String line;
                    int maxCorrida = 0;
                    while ((line = br.readLine()) != null) {
                        if (line.trim().isEmpty()) continue;
                        String[] partes = line.split(",", 2); // Solo necesitamos la primera columna
                        if (partes.length > 0) {
                            try {
                                int corrida = Integer.parseInt(partes[0].trim());
                                maxCorrida = Math.max(maxCorrida, corrida);
                            } catch (NumberFormatException e) {
                                // Ignorar líneas mal formateadas
                            }
                        }
                    }
                    corridaInicial = maxCorrida + 1;
                }
            } catch (Exception e) {
                System.err.println("Advertencia: No se pudo leer el archivo CSV existente. Empezando desde corrida 1.");
                corridaInicial = 1;
            }
        }
        
        // Abrir archivo en modo APPEND (true = agregar al final, no sobrescribir)
        try (PrintWriter csvWriter = new PrintWriter(
                new OutputStreamWriter(new FileOutputStream(RESULTADOS_GRASP_CSV, true), StandardCharsets.UTF_8))) {
            
            // Escribir encabezados solo si el archivo no existe
            if (!archivoExiste) {
                csvWriter.println("Corrida,Seed,Algoritmo,Piezas_a_posicionar,Piezas_colocadas,Area_colocada,Area_plancha,Area_desperdiciada,Ymax,Utilizacion_%,Tiempo_ms,Tiempo_s");
            }
            
            // Ejecutar algoritmo GRASP múltiples veces
            System.out.println("\n========================================");
            if (archivoExiste) {
                System.out.println("Archivo CSV existente detectado. Continuando desde corrida " + corridaInicial);
                System.out.println("Agregando " + NUM_ITERACIONES + " nuevas ejecuciones del algoritmo GRASP");
            } else {
                System.out.println("Iniciando " + NUM_ITERACIONES + " ejecuciones del algoritmo GRASP");
            }
            System.out.println("========================================");
            
            GraspResult ultimoResult = null;
            
            for (int i = 0; i < NUM_ITERACIONES; i++) {
                int corrida = corridaInicial + i;
                System.out.println("\n>>> GRASP - CORRIDA " + corrida + " (ejecución " + (i + 1) + " de " + NUM_ITERACIONES + ") <<<");
                
                // Ejecutar algoritmo GRASP y obtener estadísticas
                long seed = System.nanoTime() + corrida;
                GraspResult result = GraspPlacer.placeWithStats(sheetInset, parts, seed, GraspPlacer.DEFAULT_ALPHA);
                result.corrida = corrida;
                ultimoResult = result; // Guardar para el SVG final
                
                // Guardar resultado en CSV
                escribirResultadoGraspCSV(csvWriter, result);
                csvWriter.flush(); // Asegurar que se escriba inmediatamente
                
                System.out.println(">>> Resultado guardado en CSV: Corrida " + corrida + 
                                 " | Colocadas: " + result.placedCount + 
                                 " | Area: " + DF3.format(result.areaPlaced) + 
                                 " | Ymax: " + DF3.format(result.yMax) + 
                                 " | Utilización: " + DF3.format(result.utilizationPct) + "%");
            }
            
            System.out.println("\n========================================");
            System.out.println("Todas las corridas GRASP completadas");
            System.out.println("Resultados guardados en: " + RESULTADOS_GRASP_CSV);
            if (archivoExiste) {
                System.out.println("(Datos agregados al final del archivo existente)");
            }
            System.out.println("========================================");
            
            // Exportar SVG de la última corrida (opcional, para visualización)
            if (ultimoResult != null) {
                String svg = SvgWriter.toSVG(sheet, ultimoResult.solution);
                try (Writer w = new OutputStreamWriter(new FileOutputStream(outputSvg), StandardCharsets.UTF_8)) {
                    w.write(svg);
                }
                System.err.println("SVG generado en: " + outputSvg);
            }
        }
    }
    
    private static void escribirResultadoCSV(PrintWriter csvWriter, HybridResult result) {
        csvWriter.print(result.corrida);
        csvWriter.print(",");
        csvWriter.print(result.seed);
        csvWriter.print(",");
        csvWriter.print(result.algo != null ? result.algo : "UNKNOWN");
        csvWriter.print(",");
        csvWriter.print(result.toPlace);
        csvWriter.print(",");
        csvWriter.print(result.placedCount);
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.areaPlaced));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.areaSheet));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.areaWaste));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.yMax));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.utilizationPct));
        csvWriter.print(",");
        csvWriter.print(result.elapsedMs);
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.elapsedMs / 1000.0));
        csvWriter.println();
    }
    
    private static void escribirResultadoGraspCSV(PrintWriter csvWriter, GraspResult result) {
        csvWriter.print(result.corrida);
        csvWriter.print(",");
        csvWriter.print(result.seed);
        csvWriter.print(",");
        csvWriter.print(result.algo != null ? result.algo : "GRASP");
        csvWriter.print(",");
        csvWriter.print(result.toPlace);
        csvWriter.print(",");
        csvWriter.print(result.placedCount);
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.areaPlaced));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.areaSheet));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.areaWaste));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.yMax));
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.utilizationPct));
        csvWriter.print(",");
        csvWriter.print(result.elapsedMs);
        csvWriter.print(",");
        csvWriter.print(DF3.format(result.elapsedMs / 1000.0));
        csvWriter.println();
    }
}
