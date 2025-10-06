package org.tesis.nesting;

import org.locationtech.jts.geom.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.List;


public class TesisNesting {

    public static void main(String[] args) throws Exception {
        //Acá se puede cambiar el nombre de el svg resultante.
        String inputCsv  = "input/datosVariados.csv";
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

        // 4) Seleccionar que algoritmo se desea usar
        //Cuidado, el genético demora casi una hora.
        
        //List<PlacedPart> placed = BottomLeftPlacer.place(sheetInset, parts);
        //List<PlacedPart> placed = GraspPlacer.place(sheetInset, parts, System.nanoTime());
        //List<PlacedPart> placed = GeneticPlacer.place(sheetInset, parts);
        //List<PlacedPart> placed = MemeticPlacer.place(sheetInset, parts);
        List<PlacedPart> placed = HybridGraspMemeticPlacer.place(sheetInset, parts);
        

// 5) Exportar SVG
        String svg = SvgWriter.toSVG(sheet, placed);
        try (Writer w = new OutputStreamWriter(new FileOutputStream(outputSvg), StandardCharsets.UTF_8)) {
            w.write(svg);
        }
        System.out.println(svg);
        System.err.println("SVG generado en: " + outputSvg);
    }
}
