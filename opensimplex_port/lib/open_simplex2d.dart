// Ported from https://gist.github.com/KdotJPG/b1270127455a94ac5d19

import 'dart:math';

// Constants
const int primeX = 0x5205402B9270C86F;
const int primeY = 0x598CD327003817B5;
const int primeZ = 0x5BCC226E9FA0BACB;
const int primeW = 0x56CC5227E58F554B;

const int hashMultiplier = 0x53A3F72DEEC546F5;

const double root2over2 = 0.7071067811865476;
const double skew2d = 0.366025403784439;
const double unSkew2d = -0.21132486540518713;

const double root3over3 = 0.577350269189626;
const double fallbackRotate3 = 2.0 / 3.0;
const double rotate3orthogonalizer = unSkew2d;

const int nGrads2dExponent = 7;
const int nGrads3dExponent = 8;

const int nGrads2d = 1 << nGrads2dExponent;
const int nGrads3d = 1 << nGrads3dExponent;

const double normalizer2d = 0.05481866495625118;
const double normalizer3d = 0.2781926117527186;

const double rSquared2d = 2.0 / 3.0;
const double rSquared3d = 3.0 / 4.0;

class OpenSimplex2S {
  init() {
    initGradients();
  }

  // Noise Evaluators

  // 2D OpenSimplex2S/SuperSimplex noise, standard lattice orientation.
  double noise2(int seed, double x, double y) {
    // Get points for A2* lattice
    double s = skew2d * (x + y);
    double xs = x + s, ys = y + s;

    return noise2UnskewedBase(seed, xs, ys);
  }

  // 2D OpenSimplex2S/SuperSimplex noise, with Y pointing down the main diagonal.
  double noise2ImproveX(int seed, double x, double y) {
    // Skew transform and rotation baked into one.
    double xx = x * root2over2;
    double yy = y * (root2over2 * (1 + 2 * skew2d));

    return noise2UnskewedBase(seed, yy + xx, yy - xx);
  }

  // 2D OpenSimplex2S/SuperSimplex noise base.
  double noise2UnskewedBase(int seed, double xs, double ys) {
    // Get base points and offsets.
    int xsb = fastFloor(xs), ysb = fastFloor(ys);
    double xi = xs - xsb, yi = ys - ysb;

    // Prime pre-multiplication for hash.
    int xsbp = xsb * primeX, ysbp = ysb * primeY;

    // Unskew.
    double t = (xi + yi) * unSkew2d;
    double dx0 = xi + t, dy0 = yi + t;

    // First vertex.
    double a0 = rSquared2d - dx0 * dx0 - dy0 * dy0;
    double value = (a0 * a0) * (a0 * a0) * grad(seed, xsbp, ysbp, dx0, dy0);

    // Second vertex.
    double a1 = 2 * (1 + 2 * unSkew2d) * (1 / unSkew2d + 2) * t +
        (-2 * (1 + 2 * unSkew2d) * (1 + 2 * unSkew2d) + a0);
    double dx1 = dx0 - (1 + 2 * unSkew2d);
    double dy1 = dy0 - (1 + 2 * unSkew2d);
    value += (a1 * a1) *
        (a1 * a1) *
        grad(seed, xsbp + primeX, ysbp + primeY, dx1, dy1);

    // Third and fourth vertices.
    double xmyi = xi - yi;
    if (t < unSkew2d) {
      if (xi + xmyi > 1) {
        double dx2 = dx0 - (3 * unSkew2d + 2);
        double dy2 = dy0 - (3 * unSkew2d + 1);
        double a2 = rSquared2d - dx2 * dx2 - dy2 * dy2;
        if (a2 > 0) {
          value += (a2 * a2) *
              (a2 * a2) *
              grad(seed, xsbp + (primeX << 1), ysbp + primeY, dx2, dy2);
        }
      } else {
        double dx2 = dx0 - unSkew2d;
        double dy2 = dy0 - (unSkew2d + 1);
        double a2 = rSquared2d - dx2 * dx2 - dy2 * dy2;
        if (a2 > 0) {
          value +=
              (a2 * a2) * (a2 * a2) * grad(seed, xsbp, ysbp + primeY, dx2, dy2);
        }
      }

      if (yi - xmyi > 1) {
        double dx3 = dx0 - (3 * unSkew2d + 1);
        double dy3 = dy0 - (3 * unSkew2d + 2);
        double a3 = rSquared2d - dx3 * dx3 - dy3 * dy3;
        if (a3 > 0) {
          value += (a3 * a3) *
              (a3 * a3) *
              grad(seed, xsbp + primeX, ysbp + (primeY << 1), dx3, dy3);
        }
      } else {
        double dx3 = dx0 - (unSkew2d + 1);
        double dy3 = dy0 - unSkew2d;
        double a3 = rSquared2d - dx3 * dx3 - dy3 * dy3;
        if (a3 > 0) {
          value +=
              (a3 * a3) * (a3 * a3) * grad(seed, xsbp + primeX, ysbp, dx3, dy3);
        }
      }
    } else {
      if (xi + xmyi < 0) {
        double dx2 = dx0 + (1 + unSkew2d);
        double dy2 = dy0 + unSkew2d;
        double a2 = rSquared2d - dx2 * dx2 - dy2 * dy2;
        if (a2 > 0) {
          value +=
              (a2 * a2) * (a2 * a2) * grad(seed, xsbp - primeX, ysbp, dx2, dy2);
        }
      } else {
        double dx2 = dx0 - (unSkew2d + 1);
        double dy2 = dy0 - unSkew2d;
        double a2 = rSquared2d - dx2 * dx2 - dy2 * dy2;
        if (a2 > 0) {
          value +=
              (a2 * a2) * (a2 * a2) * grad(seed, xsbp + primeX, ysbp, dx2, dy2);
        }
      }

      if (yi < xmyi) {
        double dx2 = dx0 + unSkew2d;
        double dy2 = dy0 + (unSkew2d + 1);
        double a2 = rSquared2d - dx2 * dx2 - dy2 * dy2;
        if (a2 > 0) {
          value +=
              (a2 * a2) * (a2 * a2) * grad(seed, xsbp, ysbp - primeY, dx2, dy2);
        }
      } else {
        double dx2 = dx0 - unSkew2d;
        double dy2 = dy0 - (unSkew2d + 1);
        double a2 = rSquared2d - dx2 * dx2 - dy2 * dy2;
        if (a2 > 0) {
          value +=
              (a2 * a2) * (a2 * a2) * grad(seed, xsbp, ysbp + primeY, dx2, dy2);
        }
      }
    }

    return value;
  }

  double noise3ImproveXY(int seed, double x, double y, double z) {
    // Re-orient the cubic lattices without skewing, so Z points up the main lattice diagonal,
    // and the planes formed by XY are moved far out of alignment with the cube faces.
    // Orthonormal rotation. Not a skew transform.
    double xy = x + y;
    double s2 = xy * rotate3orthogonalizer;
    double zz = z * root3over3;
    double xr = x + s2 + zz;
    double yr = y + s2 + zz;
    double zr = xy * -root3over3 + zz;

    // Evaluate both lattices to form a BCC lattice.
    return noise3UnrotatedBase(seed, xr, yr, zr);
  }

  double noise3UnrotatedBase(int seed, double xr, double yr, double zr) {
    // Get base points and offsets.
    int xrb = fastFloor(xr);
    int yrb = fastFloor(yr);
    int zrb = fastFloor(zr);
    double xi = xr - xrb;
    double yi = yr - yrb;
    double zi = zr - zrb;

    // Prime pre-multiplication for hash. Also flip seed for second lattice copy.
    int xrbp = xrb * primeX;
    int yrbp = yrb * primeY;
    int zrbp = zrb * primeZ;
    int seed2 = seed ^ -0x52D547B2E96ED629;

    // -1 if positive, 0 if negative.
    int xNMask = (-0.5 - xi).toInt();
    int yNMask = (-0.5 - yi).toInt();
    int zNMask = (-0.5 - zi).toInt();

    // First vertex.
    double x0 = xi + xNMask;
    double y0 = yi + yNMask;
    double z0 = zi + zNMask;
    double a0 = rSquared3d - x0 * x0 - y0 * y0 - z0 * z0;
    double value = pow(a0, 4) *
        gradxyz(seed, xrbp + (xNMask & primeX), yrbp + (yNMask & primeY),
            zrbp + (zNMask & primeZ), x0, y0, z0);

    // Second vertex.
    double x1 = xi - 0.5;
    double y1 = yi - 0.5;
    double z1 = zi - 0.5;
    double a1 = rSquared3d - x1 * x1 - y1 * y1 - z1 * z1;
    value += pow(a1, 4) *
        gradxyz(seed2, xrbp + primeX, yrbp + primeY, zrbp + primeZ, x1, y1, z1);

    // Shortcuts for building the remaining falloffs.
    double xAFlipMask0 = ((xNMask | 1) << 1) * x1;
    double yAFlipMask0 = ((yNMask | 1) << 1) * y1;
    double zAFlipMask0 = ((zNMask | 1) << 1) * z1;
    double xAFlipMask1 = (-2 - (xNMask << 2)) * x1 - 1.0;
    double yAFlipMask1 = (-2 - (yNMask << 2)) * y1 - 1.0;
    double zAFlipMask1 = (-2 - (zNMask << 2)) * z1 - 1.0;

    bool skip5 = false;
    double a2 = xAFlipMask0 + a0;
    if (a2 > 0) {
      double x2 = x0 - (xNMask | 1);
      double y2 = y0;
      double z2 = z0;
      value += pow(a2, 4) *
          gradxyz(seed, xrbp + (~xNMask & primeX), yrbp + (yNMask & primeY),
              zrbp + (zNMask & primeZ), x2, y2, z2);
    } else {
      double a3 = yAFlipMask0 + zAFlipMask0 + a0;
      if (a3 > 0) {
        double x3 = x0;
        double y3 = y0 - (yNMask | 1);
        double z3 = z0 - (zNMask | 1);
        value += pow(a3, 4) *
            gradxyz(seed, xrbp + (xNMask & primeX), yrbp + (~yNMask & primeY),
                zrbp + (~zNMask & primeZ), x3, y3, z3);
      }

      double a4 = xAFlipMask1 + a1;
      if (a4 > 0) {
        double x4 = (xNMask | 1) + x1;
        double y4 = y1;
        double z4 = z1;
        value += pow(a4, 4) *
            gradxyz(seed2, xrbp + (xNMask & (primeX * 2)), yrbp + primeY,
                zrbp + primeZ, x4, y4, z4);
        skip5 = true;
      }
    }

    bool skip9 = false;
    double a6 = yAFlipMask0 + a0;
    if (a6 > 0) {
      double x6 = x0;
      double y6 = y0 - (yNMask | 1);
      double z6 = z0;
      value += pow(a6, 4) *
          gradxyz(seed, xrbp + (xNMask & primeX), yrbp + (~yNMask & primeY),
              zrbp + (zNMask & primeZ), x6, y6, z6);
    } else {
      double a7 = xAFlipMask0 + zAFlipMask0 + a0;
      if (a7 > 0) {
        double x7 = x0 - (xNMask | 1);
        double y7 = y0;
        double z7 = z0 - (zNMask | 1);
        value += pow(a7, 4) *
            gradxyz(seed, xrbp + (~xNMask & primeX), yrbp + (yNMask & primeY),
                zrbp + (~zNMask & primeZ), x7, y7, z7);
      }

      double a8 = yAFlipMask1 + a1;
      if (a8 > 0) {
        double x8 = x1;
        double y8 = (yNMask | 1) + y1;
        double z8 = z1;
        value += pow(a8, 4) *
            gradxyz(seed2, xrbp + primeX, yrbp + (yNMask & (primeY << 1)),
                zrbp + primeZ, x8, y8, z8);
        skip9 = true;
      }
    }

    bool skipD = false;
    double aA = zAFlipMask0 + a0;
    if (aA > 0) {
      double xA = x0;
      double yA = y0;
      double zA = z0 - (zNMask | 1);
      value += pow(aA, 4) *
          gradxyz(seed, xrbp + (xNMask & primeX), yrbp + (yNMask & primeY),
              zrbp + (~zNMask & primeZ), xA, yA, zA);
    } else {
      double aB = xAFlipMask0 + yAFlipMask0 + a0;
      if (aB > 0) {
        double xB = x0 - (xNMask | 1);
        double yB = y0 - (yNMask | 1);
        double zB = z0;
        value += pow(aB, 4) *
            gradxyz(seed, xrbp + (~xNMask & primeX), yrbp + (~yNMask & primeY),
                zrbp + (zNMask & primeZ), xB, yB, zB);
      }

      double aC = zAFlipMask1 + a1;
      if (aC > 0) {
        double xC = x1;
        double yC = y1;
        double zC = (zNMask | 1) + z1;
        value += pow(aC, 4) *
            gradxyz(seed2, xrbp + primeX, yrbp + primeY,
                zrbp + (zNMask & (primeZ << 1)), xC, yC, zC);
        skipD = true;
      }
    }

    if (!skip5) {
      double a5 = yAFlipMask1 + zAFlipMask1 + a1;
      if (a5 > 0) {
        double x5 = x1;
        double y5 = (yNMask | 1) + y1;
        double z5 = (zNMask | 1) + z1;
        value += pow(a5, 4) *
            gradxyz(seed2, xrbp + primeX, yrbp + (yNMask & (primeY << 1)),
                zrbp + (zNMask & (primeZ << 1)), x5, y5, z5);
      }
    }

    if (!skip9) {
      double a9 = xAFlipMask1 + zAFlipMask1 + a1;
      if (a9 > 0) {
        double x9 = (xNMask | 1) + x1;
        double y9 = y1;
        double z9 = (zNMask | 1) + z1;
        value += pow(a9, 4) *
            gradxyz(seed2, xrbp + (xNMask & (primeX * 2)), yrbp + primeY,
                zrbp + (zNMask & (primeZ << 1)), x9, y9, z9);
      }
    }

    if (!skipD) {
      double aD = xAFlipMask1 + yAFlipMask1 + a1;
      if (aD > 0) {
        double xD = (xNMask | 1) + x1;
        double yD = (yNMask | 1) + y1;
        double zD = z1;
        value += pow(aD, 4) *
            gradxyz(seed2, xrbp + (xNMask & (primeX << 1)),
                yrbp + (yNMask & (primeY << 1)), zrbp + primeZ, xD, yD, zD);
      }
    }

    return value;
  }

  // // Placeholder for the fastFloor function
  // static int fastFloor(double x) => x.floor();

  // // Placeholder for the grad function
  // static double grad(int seed, int xPrimed, int yPrimed, int zPrimed, double dx, double dy, double dz) {
  //   // Implementation of grad
  //   // This function needs to be defined based on your specific use case
  //   return 0.0; // Return appropriate value based on implementation
  // }

  // Utility
  double grad(int seed, int xsvp, int ysvp, double dx, double dy) {
    int hash = seed ^ xsvp ^ ysvp;
    hash *= hashMultiplier;
    hash ^= hash >> (64 - nGrads2dExponent + 1);
    int gi = hash & ((nGrads2d - 1) << 1);
    return gradients2d[gi | 0] * dx + gradients2d[gi | 1] * dy;
  }

  double gradxyz(
      int seed, int xrvp, int yrvp, int zrvp, double dx, double dy, double dz) {
    // Calculate hash
    int hash = (seed ^ xrvp) ^ (yrvp ^ zrvp);
    hash *= hashMultiplier;
    hash ^= hash >> (64 - nGrads3dExponent + 2);

    // Get gradient index
    int gi = hash & ((nGrads3d - 1) << 2);

    // Calculate dot product
    return gradients3d[gi | 0] * dx +
        gradients3d[gi | 1] * dy +
        gradients3d[gi | 2] * dz;
  }

  int fastFloor(double x) {
    int xi = x.toInt();
    return x < xi ? xi - 1 : xi;
  }

  // Lookup Tables & Gradients

  List<double> gradients2d = List.filled(nGrads2d * 2, 0.0);
  List<double> gradients3d = List.filled(nGrads3d * 4, 0.0);

  void initGradients() {
    List<double> grad2 = [
      0.38268343236509, 0.923879532511287,
      0.923879532511287, 0.38268343236509,
      0.923879532511287, -0.38268343236509,
      0.38268343236509, -0.923879532511287,
      -0.38268343236509, -0.923879532511287,
      -0.923879532511287, -0.38268343236509,
      -0.923879532511287, 0.38268343236509,
      -0.38268343236509, 0.923879532511287,
      //-------------------------------------//
      0.130526192220052, 0.99144486137381,
      0.608761429008721, 0.793353340291235,
      0.793353340291235, 0.608761429008721,
      0.99144486137381, 0.130526192220051,
      0.99144486137381, -0.130526192220051,
      0.793353340291235, -0.60876142900872,
      0.608761429008721, -0.793353340291235,
      0.130526192220052, -0.99144486137381,
      -0.130526192220052, -0.99144486137381,
      -0.608761429008721, -0.793353340291235,
      -0.793353340291235, -0.608761429008721,
      -0.99144486137381, -0.130526192220052,
      -0.99144486137381, 0.130526192220051,
      -0.793353340291235, 0.608761429008721,
      -0.608761429008721, 0.793353340291235,
      -0.130526192220052, 0.99144486137381
    ];

    for (int i = 0; i < grad2.length; i++) {
      grad2[i] = grad2[i] / normalizer2d;
    }

    for (int i = 0, j = 0; i < gradients2d.length; i++, j++) {
      if (j == grad2.length) j = 0;
      gradients2d[i] = grad2[j];
    }

    List<double> grad3 = [
      2.22474487139, 2.22474487139, -1.0, 0.0,
      2.22474487139, 2.22474487139, 1.0, 0.0,
      3.0862664687972017, 1.1721513422464978, 0.0, 0.0,
      1.1721513422464978, 3.0862664687972017, 0.0, 0.0,
      -2.22474487139, 2.22474487139, -1.0, 0.0,
      -2.22474487139, 2.22474487139, 1.0, 0.0,
      -1.1721513422464978, 3.0862664687972017, 0.0, 0.0,
      -3.0862664687972017, 1.1721513422464978, 0.0, 0.0,
      -1.0, -2.22474487139, -2.22474487139, 0.0,
      1.0, -2.22474487139, -2.22474487139, 0.0,
      0.0, -3.0862664687972017, -1.1721513422464978, 0.0,
      0.0, -1.1721513422464978, -3.0862664687972017, 0.0,
      -1.0, -2.22474487139, 2.22474487139, 0.0,
      1.0, -2.22474487139, 2.22474487139, 0.0,
      0.0, -1.1721513422464978, 3.0862664687972017, 0.0,
      0.0, -3.0862664687972017, 1.1721513422464978, 0.0,
      //-------------------------------------//
      -2.22474487139, -2.22474487139, -1.0, 0.0,
      -2.22474487139, -2.22474487139, 1.0, 0.0,
      -3.0862664687972017, -1.1721513422464978, 0.0, 0.0,
      -1.1721513422464978, -3.0862664687972017, 0.0, 0.0,
      -2.22474487139, -1.0, -2.22474487139, 0.0,
      -2.22474487139, 1.0, -2.22474487139, 0.0,
      -1.1721513422464978, 0.0, -3.0862664687972017, 0.0,
      -3.0862664687972017, 0.0, -1.1721513422464978, 0.0,
      -2.22474487139, -1.0, 2.22474487139, 0.0,
      -2.22474487139, 1.0, 2.22474487139, 0.0,
      -3.0862664687972017, 0.0, 1.1721513422464978, 0.0,
      -1.1721513422464978, 0.0, 3.0862664687972017, 0.0,
      -1.0, 2.22474487139, -2.22474487139, 0.0,
      1.0, 2.22474487139, -2.22474487139, 0.0,
      0.0, 1.1721513422464978, -3.0862664687972017, 0.0,
      0.0, 3.0862664687972017, -1.1721513422464978, 0.0,
      -1.0, 2.22474487139, 2.22474487139, 0.0,
      1.0, 2.22474487139, 2.22474487139, 0.0,
      0.0, 3.0862664687972017, 1.1721513422464978, 0.0,
      0.0, 1.1721513422464978, 3.0862664687972017, 0.0,
      2.22474487139, -2.22474487139, -1.0, 0.0,
      2.22474487139, -2.22474487139, 1.0, 0.0,
      1.1721513422464978, -3.0862664687972017, 0.0, 0.0,
      3.0862664687972017, -1.1721513422464978, 0.0, 0.0,
      2.22474487139, -1.0, -2.22474487139, 0.0,
      2.22474487139, 1.0, -2.22474487139, 0.0,
      3.0862664687972017, 0.0, -1.1721513422464978, 0.0,
      1.1721513422464978, 0.0, -3.0862664687972017, 0.0,
      2.22474487139, -1.0, 2.22474487139, 0.0,
      2.22474487139, 1.0, 2.22474487139, 0.0,
      1.1721513422464978, 0.0, 3.0862664687972017, 0.0,
      3.0862664687972017, 0.0, 1.1721513422464978, 0.0,
    ];

    for (int i = 0; i < grad3.length; i++) {
      grad3[i] = grad3[i] / normalizer3d;
    }

    for (int i = 0, j = 0; i < gradients3d.length; i++, j++) {
      if (j == grad3.length) j = 0;
      gradients3d[i] = grad3[j];
    }
  }
}
