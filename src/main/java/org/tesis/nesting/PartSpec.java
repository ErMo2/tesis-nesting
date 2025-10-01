package org.tesis.nesting;

import org.locationtech.jts.geom.Polygon;

public class PartSpec {
    String id;
    Polygon poly;
    int     qty = 1;
    int     rotationStepDeg = 1;
    boolean mirrorAllowed   = false;
    double  clearance       = 0.0;
    double  grainAngleDeg   = 0.0;
}
