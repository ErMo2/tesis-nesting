package org.tesis.nesting;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

public class CsvGeometryReader {

    // función para leer un CSV y devolver la lista de vértices
    static List<VertexRow> readCsv(String path) throws IOException {
        List<VertexRow> out = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(path), StandardCharsets.UTF_8))) {
            String header = br.readLine();
            if (header == null) throw new IOException("CSV vacío: " + path);
            String[] h = VertexRow.splitCsv(header);
            String line;
            while ((line = br.readLine()) != null) {
                if (line.trim().isEmpty()) continue;
                String[] v = VertexRow.splitCsv(line);
                out.add(VertexRow.fromCsv(h, v));
            }
        }
        return out;
    }
}
