// Constants
const int primeX = 0x5205402B9270C86F;
const int primeY = 0x598CD327003817B5;
const int hashMultiplier = 0x53A3F72DEEC546F5;

const double root2over2 = 0.7071067811865476;
const double skew2d = 0.366025403784439;
const double unSkew2d = -0.21132486540518713;

const int nGrads2dExponent = 7;
const int nGrads2d = 1 << nGrads2dExponent;

const double normalizer2d = 0.05481866495625118;

const double rSquared2d = 2.0 / 3.0;

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

// Utility
double grad(int seed, int xsvp, int ysvp, double dx, double dy) {
  int hash = seed ^ xsvp ^ ysvp;
  hash *= hashMultiplier;
  hash ^= hash >> (64 - nGrads2dExponent + 1);
  int gi = hash & ((nGrads2d - 1) << 1);
  return gradients2d[gi | 0] * dx + gradients2d[gi | 1] * dy;
}

int fastFloor(double x) {
  int xi = x.toInt();
  return x < xi ? xi - 1 : xi;
}

// Lookup Tables & Gradients
List<double> gradients2d = List.filled(nGrads2d * 2, 0.0);

void init() {
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
}
