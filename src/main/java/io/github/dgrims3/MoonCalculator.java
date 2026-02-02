package io.github.dgrims3;

import java.time.LocalDate;
import java.time.LocalDateTime;
import java.time.ZonedDateTime;
import java.util.ArrayList;
import java.util.List;

public abstract class MoonCalculator extends AbstractAstronomicalCalculator {
    // moonLat β, moonLong λ, earthMoonDistance ∆

    int COUNT = 60;

    enum position {ASCENSION, DECLINATION}

    enum time {Rise, Transit, Set}

    enum coords {LAT, LNG, DISTANCE}

    enum Term {
        /**
         * Mean Elongation of the Moon.
         */
        TERM_D,
        /**
         * Sun's Mean Anomaly
         */
        TERM_M,
        /**
         * Moon;s Mean Anomaly
         */
        TERM_MPR,
        /**
         * Moon's Mean Argument of Latitude (mean distance of the moon from it's ascending node).
         */
        TERM_F,
        /**
         * Coefficient of the sine of the Argument
         */
        TERM_LB,
        /**
         * Coefficient of the cosine of the Argument
         */
        TERM_R
    }

    //  Moon's Periodic Terms for Longitude and Distance.
    double[][] ML_TERMS =
            {
                    {0, 0, 1, 0, 6288774, -20905355},
                    {2, 0, -1, 0, 1274027, -3699111},
                    {2, 0, 0, 0, 658314, -2955968},
                    {0, 0, 2, 0, 213618, -569925},
                    {0, 1, 0, 0, -185116, 48888},
                    {0, 0, 0, 2, -114332, -3149},
                    {2, 0, -2, 0, 58793, 246158},
                    {2, -1, -1, 0, 57066, -152138},
                    {2, 0, 1, 0, 53322, -170733},
                    {2, -1, 0, 0, 45758, -204586},
                    {0, 1, -1, 0, -40923, -129620},
                    {1, 0, 0, 0, -34720, 108743},
                    {0, 1, 1, 0, -30383, 104755},
                    {2, 0, 0, -2, 15327, 10321},
                    {0, 0, 1, 2, -12528, 0},
                    {0, 0, 1, -2, 10980, 79661},
                    {4, 0, -1, 0, 10675, -34782},
                    {0, 0, 3, 0, 10034, -23210},
                    {4, 0, -2, 0, 8548, -21636},
                    {2, 1, -1, 0, -7888, 24208},
                    {2, 1, 0, 0, -6766, 30824},
                    {1, 0, -1, 0, -5163, -8379},
                    {1, 1, 0, 0, 4987, -16675},
                    {2, -1, 1, 0, 4036, -12831},
                    {2, 0, 2, 0, 3994, -10445},
                    {4, 0, 0, 0, 3861, -11650},
                    {2, 0, -3, 0, 3665, 14403},
                    {0, 1, -2, 0, -2689, -7003},
                    {2, 0, -1, 2, -2602, 0},
                    {2, -1, -2, 0, 2390, 10056},
                    {1, 0, 1, 0, -2348, 6322},
                    {2, -2, 0, 0, 2236, -9884},
                    {0, 1, 2, 0, -2120, 5751},
                    {0, 2, 0, 0, -2069, 0},
                    {2, -2, -1, 0, 2048, -4950},
                    {2, 0, 1, -2, -1773, 4130},
                    {2, 0, 0, 2, -1595, 0},
                    {4, -1, -1, 0, 1215, -3958},
                    {0, 0, 2, 2, -1110, 0},
                    {3, 0, -1, 0, -892, 3258},
                    {2, 1, 1, 0, -810, 2616},
                    {4, -1, -2, 0, 759, -1897},
                    {0, 2, -1, 0, -713, -2117},
                    {2, 2, -1, 0, -700, 2354},
                    {2, 1, -2, 0, 691, 0},
                    {2, -1, 0, -2, 596, 0},
                    {4, 0, 1, 0, 549, -1423},
                    {0, 0, 4, 0, 537, -1117},
                    {4, -1, 0, 0, 520, -1571},
                    {1, 0, -2, 0, -487, -1739},
                    {2, 1, 0, -2, -399, 0},
                    {0, 0, 2, -2, -381, -4421},
                    {1, 1, 1, 0, 351, 0},
                    {3, 0, -2, 0, -340, 0},
                    {4, 0, -3, 0, 330, 0},
                    {2, -1, 2, 0, 327, 0},
                    {0, 2, 1, 0, -323, 1165},
                    {1, 1, -1, 0, 299, 0},
                    {2, 0, 3, 0, 294, 0},
                    {2, 0, -1, -2, 0, 8752}
            };

    //  Moon's Periodic Terms for Latitude
    double[][] MB_TERMS =
            {
                    {0, 0, 0, 1, 5128122, 0},
                    {0, 0, 1, 1, 280602, 0},
                    {0, 0, 1, -1, 277693, 0},
                    {2, 0, 0, -1, 173237, 0},
                    {2, 0, -1, 1, 55413, 0},
                    {2, 0, -1, -1, 46271, 0},
                    {2, 0, 0, 1, 32573, 0},
                    {0, 0, 2, 1, 17198, 0},
                    {2, 0, 1, -1, 9266, 0},
                    {0, 0, 2, -1, 8822, 0},
                    {2, -1, 0, -1, 8216, 0},
                    {2, 0, -2, -1, 4324, 0},
                    {2, 0, 1, 1, 4200, 0},
                    {2, 1, 0, -1, -3359, 0},
                    {2, -1, -1, 1, 2463, 0},
                    {2, -1, 0, 1, 2211, 0},
                    {2, -1, -1, -1, 2065, 0},
                    {0, 1, -1, -1, -1870, 0},
                    {4, 0, -1, -1, 1828, 0},
                    {0, 1, 0, 1, -1794, 0},
                    {0, 0, 0, 3, -1749, 0},
                    {0, 1, -1, 1, -1565, 0},
                    {1, 0, 0, 1, -1491, 0},
                    {0, 1, 1, 1, -1475, 0},
                    {0, 1, 1, -1, -1410, 0},
                    {0, 1, 0, -1, -1344, 0},
                    {1, 0, 0, -1, -1335, 0},
                    {0, 0, 3, 1, 1107, 0},
                    {4, 0, 0, -1, 1021, 0},
                    {4, 0, -1, 1, 833, 0},
                    {0, 0, 1, -3, 777, 0},
                    {4, 0, -2, 1, 671, 0},
                    {2, 0, 0, -3, 607, 0},
                    {2, 0, 2, -1, 596, 0},
                    {2, -1, 1, -1, 491, 0},
                    {2, 0, -2, 1, -451, 0},
                    {0, 0, 3, -1, 439, 0},
                    {2, 0, 2, 1, 422, 0},
                    {2, 0, -3, -1, 421, 0},
                    {2, 1, -1, 1, -366, 0},
                    {2, 1, 0, 1, -351, 0},
                    {4, 0, 0, 1, 331, 0},
                    {2, -1, 1, 1, 315, 0},
                    {2, -2, 0, -1, 302, 0},
                    {0, 0, 1, 3, -283, 0},
                    {2, 1, 1, -1, -229, 0},
                    {1, 1, 0, -1, 223, 0},
                    {1, 1, 0, 1, 223, 0},
                    {0, 1, -2, -1, -220, 0},
                    {2, 1, -1, -1, -220, 0},
                    {1, 0, 1, 1, -185, 0},
                    {2, -1, -2, -1, 181, 0},
                    {0, 1, 2, 1, -177, 0},
                    {4, 0, -2, -1, 176, 0},
                    {4, -1, -1, -1, 166, 0},
                    {1, 0, 1, -1, -164, 0},
                    {4, 0, 1, -1, 132, 0},
                    {1, 0, -1, -1, -119, 0},
                    {4, -1, 0, -1, 115, 0},
                    {2, -2, 0, 1, 107, 0}
            };

    /**
     * Evaluates a third-order (cubic) polynomial using Horner's method for efficient computation.
     * The polynomial is of the form: f(x) = ax³ + bx² + cx + d
     * This method is commonly used to calculate astronomical parameters that vary with time,
     * where the time variable is measured in Julian centuries.
     *
     * @param a The coefficient of the cubic term (x³)
     * @param b The coefficient of the quadratic term (x²)
     * @param c The coefficient of the linear term (x)
     * @param d The constant term
     * @param jce Julian centuries from the J2000.0 epoch (the independent variable)
     * @return The evaluated polynomial value
     */
    public double third_order_polynomial(double a, double b, double c, double d, double jce) {
        return ((a * jce + b) * jce + c) * jce + d;
    }

    /**
     * Evaluates a fourth-order (quartic) polynomial using Horner's method for efficient computation.
     * The polynomial is of the form: f(x) = ax⁴ + bx³ + cx² + dx + e
     * This method is commonly used to calculate high-precision astronomical parameters that vary
     * with time, where the time variable is measured in Julian centuries. The fourth-order terms
     * account for long-term variations in lunar and planetary orbital elements.

     * @param a The coefficient of the quartic term (x⁴)
     * @param b The coefficient of the cubic term (x³)
     * @param c The coefficient of the quadratic term (x²)
     * @param d The coefficient of the linear term (x)
     * @param e The constant term
     * @param jce Julian centuries from the J2000.0 epoch (the independent variable)
     * @return The evaluated polynomial value
     */
    public double fourth_order_polynomial(double a, double b, double c, double d, double e, double jce) {
        return (((a * jce + b) * jce + c) * jce + d) * jce + e;
    }

    /**
     * Calculates the Moon's mean longitude (L'), which is the Moon's average ecliptic longitude
     * as if it moved in a circular orbit. This is one of the fundamental arguments used in
     * lunar position calculations. The mean longitude increases by approximately 13.18 degrees
     * per day as the Moon orbits Earth.

     * The calculation uses a fourth-order polynomial to account for long-term variations
     * in the Moon's orbital motion, including perturbations from the Sun and other planets.
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The Moon's mean longitude in degrees, normalized to the range [0, 360)
     */
    public double moon_mean_longitude_L_PRIME(double jce) {
        return limit_degrees(fourth_order_polynomial(
                -1.0 / 65194000, 1.0 / 538841, -0.0015786, 481267.88123421, 218.3164477, jce));
    }

    /**
     * Calculates the Moon's mean elongation (D), which is the angular separation between
     * the Moon and the Sun as seen from Earth, measured along the ecliptic.
     * When D = 0°, it's a new moon; when D = 180°, it's a full moon.

     * This parameter is crucial for determining lunar phases and is one of the fundamental
     * arguments in lunar theory. It increases by approximately 12.19 degrees per day.

     * @param jce Julian centuries from the J2000.0 epoch
     * @return The Moon's mean elongation in degrees, normalized to the range [0, 360)
     */
    public double moon_mean_elongation_D(double jce) {
        return limit_degrees(fourth_order_polynomial(
                -1.0 / 113065000, 1.0 / 545868, -0.0018819, 445267.1114034, 297.8501921, jce));
    }

    /**
     * Calculates the Sun's mean anomaly (M), which is the angular distance of the Sun from
     * its perigee (point of closest approach to Earth) as measured along its elliptical orbit.
     * When M = 0°, the Sun is at perigee; when M = 180°, it's at apogee (farthest point).

     * The Sun's mean anomaly affects the Moon's position through gravitational perturbations.
     * It increases by approximately 0.99 degrees per day, completing one cycle per year.
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The Sun's mean anomaly in degrees, normalized to the range [0, 360)
     */
    public double sun_mean_anomaly_M(double jce) {
        return limit_degrees(third_order_polynomial(
                1.0 / 24490000, -0.0001536, 35999.0502909, 357.5291092, jce));
    }

    /**
     * Calculates the Moon's mean anomaly (M'), which is the angular distance of the Moon from
     * its perigee (point of closest approach to Earth) as measured along its elliptical orbit.
     * When M' = 0°, the Moon is at perigee; when M' = 180°, it's at apogee (farthest point).

     * The Moon's mean anomaly is critical for calculating the Moon's varying speed and distance
     * from Earth. It increases by approximately 13.06 degrees per day, completing one full
     * cycle approximately every 27.55 days (the anomalistic month).
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The Moon's mean anomaly in degrees, normalized to the range [0, 360)
     */
    public double moon_mean_anomaly_M_PRIME(double jce) {
        return limit_degrees(fourth_order_polynomial(
                -1.0 / 14712000, 1.0 / 69699, 0.0087414, 477198.8675055, 134.9633964, jce));
    }

    /**
     * Calculates the Moon's argument of latitude (F), which is the mean angular distance
     * of the Moon from its ascending node along its orbit. The ascending node is the point
     * where the Moon's orbit crosses the ecliptic plane from south to north.

     * This parameter is essential for calculating the Moon's ecliptic latitude and predicting
     * lunar and solar eclipses. When F = 0° or 180°, the Moon crosses the ecliptic plane,
     * and eclipses are possible if other conditions align. It increases by approximately
     * 13.23 degrees per day, completing one cycle every 27.21 days (the draconic month).
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The Moon's argument of latitude in degrees, normalized to the range [0, 360)
     */
    public double moon_latitude_argument_F(double jce) {
        return limit_degrees(fourth_order_polynomial(
                1.0 / 863310000, -1.0 / 3526000, -0.0036539, 483202.0175233, 93.2720950, jce));
    }

    /**
     * Calculates the first additive term (A1) used in lunar longitude and latitude corrections.
     * This term accounts for the action of Venus on the Moon's position, representing one of
     * the perturbative effects from other planets in the solar system.

     * Although Venus is much smaller than the Sun, its gravitational influence still creates
     * measurable variations in the Moon's orbital elements over time.
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The A1 correction angle in degrees, normalized to the range [0, 360)
     */
    public double summationAdditiveA1(double jce) {
        return limit_degrees(119.75 + 131.859 * jce);
    }

    /**
     * Calculates the second additive term (A2) used in lunar longitude corrections.
     * This term accounts for the action of Jupiter on the Moon's position, representing
     * another planetary perturbation effect in lunar orbital calculations.

     * Jupiter, being the most massive planet in the solar system, exerts a significant
     * gravitational influence that affects the Moon's orbit, though much smaller than
     * the Sun's effect.
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The A2 correction angle in degrees, normalized to the range [0, 360)
     */
    public double summationAdditiveA2(double jce) {
        return limit_degrees(53.09 + 479264.290 * jce);
    }

    /**
     * Calculates the third additive term (A3) used in lunar latitude corrections.
     * This term accounts for Earth's orbital flattening and the oblateness of the terrestrial
     * orbit, representing additional perturbations in the Moon's position calculations.

     * This correction term helps achieve higher accuracy in determining the Moon's ecliptic
     * latitude by accounting for subtle variations in Earth's shape and orbit.
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The A3 correction angle in degrees, normalized to the range [0, 360)
     */
    public double summationAdditiveA3(double jce) {
        return limit_degrees(313.45 + 481266.484 * jce);
    }

    /**
     * Calculates additional correction terms to be added to the Moon's ecliptic latitude.
     * These corrections account for various perturbative effects from the Sun, planets,
     * and Earth's shape on the Moon's position perpendicular to the ecliptic plane.

     * The method combines six periodic terms with different amplitudes and periods,
     * each representing a specific gravitational or orbital perturbation. These small
     * corrections are essential for achieving arc-second precision in lunar positioning.
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The sum of latitude correction terms in units of 10^-6 degrees (microarcseconds)
     */
    public double additiveToMoonLat(double jce) {
        return ((-2235 * Math.sin(moon_mean_longitude_L_PRIME(jce))) +
                (382 * Math.sin(degToRad(summationAdditiveA3(jce)))) +
                (175 * Math.sin(degToRad(summationAdditiveA1(jce) - moon_latitude_argument_F(jce)))) +
                (175 * Math.sin(degToRad(summationAdditiveA1(jce) + moon_latitude_argument_F(jce)))) +
                (127 * Math.sin(degToRad(moon_mean_longitude_L_PRIME(jce) - moon_mean_anomaly_M_PRIME(jce)))) -
                (115 * Math.sin(degToRad(moon_mean_longitude_L_PRIME(jce) + moon_mean_anomaly_M_PRIME(jce)))));
    }

    /**
     * Calculates additional correction terms to be added to the Moon's ecliptic longitude.
     * These corrections account for various perturbative effects from Venus, Jupiter, and
     * the Moon's orbital characteristics that cause deviations from the mean longitude.

     * The method combines three periodic terms representing different sources of perturbation.
     * These corrections help refine the Moon's position along the ecliptic to achieve
     * higher accuracy in longitude calculations.
     *
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The sum of longitude correction terms in units of 10^-6 degrees (microarcseconds)
     */
    public double additiveToMoonLng(double jce) {
        return ((3958 * Math.sin(degToRad(summationAdditiveA1(jce)))) +
                (1962 * Math.sin(degToRad(moon_mean_longitude_L_PRIME(jce) - moon_latitude_argument_F(jce)))) +
                (318 * Math.sin(degToRad(summationAdditiveA2(jce)))));
    }

    /**
     * Calculates the summation of periodic terms used in lunar position theory.
     * This method implements the core of Chapront's ELP-2000/82 lunar theory, which uses
     * periodic terms (combinations of fundamental arguments) to model the complex gravitational
     * interactions affecting the Moon's orbit.

     * Each term in the summation represents a specific periodic perturbation caused by
     * gravitational forces from the Sun, planets, and orbital characteristics. The terms
     * are combined with sine and cosine functions weighted by an eccentricity correction factor.
     *
     * @param d Moon's mean elongation (angular distance from the Sun) in degrees
     * @param m Sun's mean anomaly in degrees
     * @param m_prime Moon's mean anomaly in degrees
     * @param f Moon's argument of latitude in degrees
     * @param jce Julian centuries from the J2000.0 epoch
     * @param terms Array of periodic term coefficients containing D, M, M', F multipliers
     *              and amplitudes for longitude/latitude and distance calculations
     * @return Array containing [sine summation for longitude/latitude, cosine summation for distance]
     */
    public double[] moon_periodic_term_summation(double d, double m, double m_prime, double f, double jce, double[][] terms) {
        double e, e_multi, trig_arg, sin_sum, cos_sum;
        e = 1.0 - jce * (0.002516 + jce * 0.0000074);
        sin_sum = 0;
        cos_sum = 0;
        for (int i = 0; i < COUNT; i++) {
            e_multi = Math.pow(e, Math.abs(terms[i][Term.TERM_M.ordinal()]));
            trig_arg = degToRad(terms[i][Term.TERM_D.ordinal()] * d + terms[i][Term.TERM_M.ordinal()] * m + terms[i][Term.TERM_MPR.ordinal()] * m_prime + terms[i][Term.TERM_F.ordinal()] * f);
            sin_sum += e_multi * terms[i][Term.TERM_LB.ordinal()] * Math.sin(trig_arg);
            cos_sum += e_multi * terms[i][Term.TERM_R.ordinal()] * Math.cos(trig_arg);
        }
        return new double[]{sin_sum, cos_sum};
    }

    /**
     * Converts the Moon's distance from its mean value plus corrections into kilometers.
     * The Moon's orbit is elliptical, so its distance from Earth varies between approximately
     * 356,500 km (perigee) and 406,700 km (apogee), with a mean distance of about 385,000 km.

     * This method takes the correction term calculated from periodic perturbations and
     * adds it to the mean Earth-Moon distance to get the actual distance at a given time.
     *
     * @param dist The distance correction term in meters from the periodic term summation
     * @return The actual Earth-Moon distance in kilometers
     */
    public double moon_earth_distance_in_km(double dist) {
        return 385000.56 + dist / 1000;
    }

    /**
     * Calculates the Moon's true ecliptic longitude by combining the mean longitude with
     * perturbation corrections. The ecliptic longitude measures the Moon's position along
     * the ecliptic plane (the plane of Earth's orbit around the Sun), starting from the
     * vernal equinox point.

     * This method adds correction terms (measured in microarcseconds) to the mean longitude
     * to account for gravitational perturbations from the Sun, planets, and orbital mechanics.
     *
     * @param longitude The longitude correction from periodic term summation (in microarcseconds)
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The Moon's true ecliptic longitude in degrees
     */
    public double moon_longitude_coordinates_in_degrees(double longitude, double jce) {
        return moon_mean_longitude_L_PRIME(jce) + (longitude + additiveToMoonLng(jce)) / 1000000;
    }

    /**
     * Calculates the Moon's ecliptic latitude by combining periodic corrections with
     * additional perturbation terms. The ecliptic latitude measures the Moon's angular
     * distance perpendicular to the ecliptic plane (above or below it).

     * The Moon's orbit is inclined about 5.1 degrees to the ecliptic, causing its latitude
     * to oscillate between approximately +5.1° and -5.1° over the course of one orbital cycle.
     * This inclination is why lunar eclipses don't occur every month.
     *
     * @param latitude The latitude correction from periodic term summation (in microarcseconds)
     * @param jce Julian centuries from the J2000.0 epoch
     * @return The Moon's ecliptic latitude in degrees
     */
    public double moon_latitude_coordinates_in_degrees(double latitude, double jce) {
        return (latitude + additiveToMoonLat(jce)) / 1000000;
    }

    /**
     * Calculates the Moon's equatorial horizontal parallax, which is the angular difference
     * in the Moon's apparent position when viewed from Earth's center versus from a point
     * on Earth's equator. This parallax is a measure of how the Moon's position shifts
     * due to the observer's location on Earth's surface.

     * Parallax is inversely proportional to distance, so it's larger when the Moon is closer
     * (at perigee, about 61 arcminutes) and smaller when farther (at apogee, about 54 arcminutes).
     * This value is crucial for accurate topocentric (observer-specific) position calculations.
     *
     * @param distance The Earth-Moon distance in kilometers
     * @return The equatorial horizontal parallax in radians
     */
    public double equatorialHorizontalParallax(double distance) {
        return Math.asin(6378.14 / distance);
    }

    /**
     * Calculates the Moon's geocentric (Earth-centered) ecliptic coordinates at a given time.
     * This method computes the Moon's position in the ecliptic coordinate system, which uses
     * the plane of Earth's orbit around the Sun as its reference plane.

     * The calculation uses Chapront's ELP-2000/82 lunar theory with periodic term summations
     * for both latitude (MB_TERMS) and longitude/distance (ML_TERMS). The fundamental arguments
     * (D, M, M', F) are combined with term coefficients to produce high-accuracy positions.
     *
     * @param jd Julian day number representing the date and time of calculation
     * @return Array containing [ecliptic latitude in degrees, ecliptic longitude in degrees,
     *         Earth-Moon distance in kilometers]
     */
    // page 340/351
    public double[] getMoonLatLngDist(double jd) {
        double jce = calcTimeJulianCent(jd);
        double lat = moon_periodic_term_summation(
                moon_mean_elongation_D(jce),
                sun_mean_anomaly_M(jce),
                moon_mean_anomaly_M_PRIME(jce),
                moon_latitude_argument_F(jce),
                jce,
                MB_TERMS
        )[0];
        double[] lngDist = moon_periodic_term_summation(
                moon_mean_elongation_D(jce),
                sun_mean_anomaly_M(jce),
                moon_mean_anomaly_M_PRIME(jce),
                moon_latitude_argument_F(jce),
                jce,
                ML_TERMS
        );
        return new double[]{moon_latitude_coordinates_in_degrees(lat, jce), moon_longitude_coordinates_in_degrees(lngDist[0], jce), moon_earth_distance_in_km(lngDist[1])};
    }

    /**
     * Converts the Moon's position from ecliptic coordinates to equatorial coordinates.
     * Equatorial coordinates (right ascension and declination) are based on Earth's equator
     * and celestial equator, making them more useful for determining where the Moon appears
     * in the sky from a specific location on Earth.

     * Right ascension (α) is measured eastward along the celestial equator from the vernal equinox,
     * similar to longitude on Earth. Declination (δ) is the angular distance north or south
     * of the celestial equator, similar to latitude. These coordinates rotate with Earth,
     * making them ideal for calculating rise, transit, and set times.
     *
     * @param jd Julian day number representing the date and time of calculation
     * @return Array containing [right ascension in degrees, declination in degrees,
     *         Earth-Moon distance in kilometers]
     */
    public double[] lunarAscensionDeclinationDistance(double jd) {
        double jce = calcTimeJulianCent(jd);
        double[] moonLatLngDist = getMoonLatLngDist(jd);
        double moonDeclination = getDeclination(moonLatLngDist[coords.LAT.ordinal()], moonLatLngDist[coords.LNG.ordinal()], jce);
        double moonRightAscension = getRightAscension(moonLatLngDist[coords.LAT.ordinal()], moonLatLngDist[coords.LNG.ordinal()], jce);
        return new double[]{moonRightAscension, moonDeclination, moonLatLngDist[coords.DISTANCE.ordinal()]};
    }

    /**
     * Calculates the Moon's geographic latitude at a specific instant during the day.
     * This is the latitude of the Moon's sub-point on Earth (where the Moon appears
     * directly overhead). The geographic latitude equals the Moon's declination.
     *
     * @param jd Julian day number at 0h UT for the date
     * @param fractionOfDay Fraction of the day (0.0 = midnight, 0.5 = noon, 1.0 = next midnight)
     * @return The Moon's geographic latitude in degrees (positive north, negative south)
     */
    public double getMoonLatAtInstant(double jd, double fractionOfDay) {
        double[] ascDecDistance = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        return ascDecDistance[position.DECLINATION.ordinal()];
    }

    /**
     * Calculates the Moon's geographic longitude at a specific instant during the day.
     * This is the longitude of the Moon's sub-point on Earth (where the Moon appears
     * directly overhead). The longitude is calculated by comparing the Moon's right
     * ascension to Greenwich Mean Sidereal Time.
     *
     * @param jd Julian day number at 0h UT for the date
     * @param fractionOfDay Fraction of the day (0.0 = midnight, 0.5 = noon, 1.0 = next midnight)
     * @return The Moon's geographic longitude in degrees (positive east, negative west)
     */
    public double getMoonLongAtInstant(double jd, double fractionOfDay) {
        double jce = calcTimeJulianCent(jd);
        double gmst = greenwichMeanSiderealTime(jce);
        double gmstAtInstant = siderealTimeAtInstantAtGreenwichInDegrees(gmst, fractionOfDay);
        double[] ascDecDistance = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        return -1 * (gmstAtInstant - ascDecDistance[position.ASCENSION.ordinal()]);
    }

    /**
     * Calculates the Moon's sub-point (the point on Earth directly beneath the Moon) at a
     * specific instant during the day. The sub-point is expressed in geographic coordinates
     * (latitude and longitude), showing where the Moon would appear directly overhead.

     * This calculation converts the Moon's equatorial coordinates to geographic coordinates
     * by accounting for Earth's rotation. The longitude is determined by comparing the Moon's
     * right ascension to Greenwich Mean Sidereal Time, effectively projecting the Moon's
     * celestial position onto Earth's surface.
     *
     * @param jd Julian day number at 0h UT for the date
     * @param fractionOfDay Fraction of the day (0.0 = midnight, 0.5 = noon, 1.0 = next midnight)
     * @return Array containing [geographic latitude in degrees, geographic longitude in degrees,
     *         Earth-Moon distance in kilometers]
     */
    public double[] getMoonLatLongDistanceAtInstant(double jd, double fractionOfDay) {
        double lat = getMoonLatAtInstant(jd, fractionOfDay);
        double lng = getMoonLongAtInstant(jd, fractionOfDay);
        double[] ascDecDistance = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        double distance = ascDecDistance[coords.DISTANCE.ordinal()];
        return new double[]{lat, lng, distance};
    }

    /**
     * Calculates the Moon's horizontal coordinates (azimuth and altitude) as seen from a
     * specific location on Earth at a given time. These coordinates describe where to look
     * in the sky to observe the Moon.

     * Azimuth is the compass direction measured clockwise from north (0° = north, 90° = east,
     * 180° = south, 270° = west). Altitude (or elevation) is the angle above the horizon
     * (0° = horizon, 90° = directly overhead at zenith, negative values = below horizon).

     * The calculation converts the Moon's equatorial coordinates to horizontal coordinates
     * using the observer's geographic location and the local sidereal time, which accounts
     * for Earth's rotation.
     *
     * @param localDateTime The local date and time of observation
     * @param offset Time zone offset from UTC in hours (e.g., -5 for EST, +1 for CET)
     * @param observersLatitude Observer's latitude in degrees (positive north, negative south)
     * @param observersLongitude Observer's longitude in degrees (positive east, negative west)
     * @return Array containing [azimuth in degrees (0-360), altitude in degrees (-90 to +90)]
     */
    public double[] getAzimuthAndAltitudeForMoonAtInstant(LocalDateTime localDateTime, int offset, double observersLatitude, double observersLongitude) {
        localDateTime = localDateTime.minusHours(offset);
        double fractionOfDay = localDateTimeToFractionOfDay(localDateTime);
        double jd = getJDFromCalenderDate(localDateTime.getYear(), localDateTime.getMonthValue(), localDateTime.getDayOfMonth());
        double sidereal = greenwichApparentSiderealTime(calcTimeJulianCent(jd));
        double siderealAtInstant = siderealTimeAtInstantAtGreenwichInDegrees(sidereal, fractionOfDay);
        double[] ascDec = lunarAscensionDeclinationDistance(jd + fractionOfDay);
        double localHourAngle = localHourAngle(siderealAtInstant,  -1 * observersLongitude, ascDec[position.ASCENSION.ordinal()]);
        double azimuth = azimuth(localHourAngle, observersLatitude, ascDec[position.DECLINATION.ordinal()]);
        azimuth += 180;
        if (azimuth == 360) {
            azimuth = 0;
        }
        double altitude = localAltitude(localHourAngle, observersLatitude, ascDec[position.DECLINATION.ordinal()]);
        return new double[]{azimuth, altitude};
    }

    /**
     * Calculates the illuminated fraction of the Moon's visible disk as a percentage.
     * This represents how much of the Moon appears lit by the Sun as seen from Earth,
     * which determines the lunar phase.

     * The calculation determines the Moon's phase angle (the angle between the Sun, Moon,
     * and Earth) and uses it to compute the illuminated fraction. Values range from 0%
     * (new moon, completely dark) to 100% (full moon, completely illuminated). Quarter
     * moons are approximately 50% illuminated.

     * The phase angle is refined using several periodic correction terms that account for
     * the elliptical nature of both the Moon's and Earth's orbits, providing accurate
     * illumination values throughout the lunar cycle.
     *
     * @param jd Julian day number representing the date and time of calculation
     * @return The percentage of the Moon's disk that is illuminated (0.0 to 100.0)
     */
    public double moonIlluminatedFractionOfDisk (double jd) {
        double jce = calcTimeJulianCent(jd);
        double D = moon_mean_elongation_D(jce);
        double M = sun_mean_anomaly_M(jce);
        double MPR = moon_mean_anomaly_M_PRIME(jce);

        double phaseAngle = 180 - D -
                (6.289 * Math.sin(MPR)) +
                (2.1 * Math.sin(M)) -
                (1.274 * Math.sin((2 * D)-M)) -
                (.658 + Math.sin(2 * D)) -
                (.214 * Math.sin(2 * MPR)) -
                (.11 * Math.sin(D));
        return ((1 + Math.cos(degToRad(phaseAngle))) / 2) * 100;
    }

    /**
     * Calculates approximate times for the Moon's rising, transit (highest point), and setting
     * at a given location. These are initial estimates that should be iteratively refined
     * for higher accuracy, as the Moon's position changes significantly during a day.

     * Rise time is when the Moon crosses the horizon ascending (altitude = 0°), transit is
     * when it crosses the observer's meridian (reaching maximum altitude), and set time is
     * when it crosses the horizon descending. The calculation accounts for atmospheric
     * refraction (approximately 34 arcminutes) and horizontal parallax.

     * The standard altitude used (-0.583°) accounts for the Moon's semi-diameter and
     * atmospheric refraction, so the calculated rise/set times represent when the Moon's
     * upper limb first becomes visible or disappears.
     *
     * @param jd Julian day number at 0h UT for the date
     * @param latlng Array containing [observer's latitude in degrees, observer's longitude in degrees]
     * @param lunarAscDecDist Array containing [right ascension in degrees, declination in degrees,
     *                        Earth-Moon distance in kilometers]
     * @return Array containing [rise time, transit time, set time] as fractions of a day (0.0 to 1.0)
     */
    public double[] risingTransitSettingApproximate(double jd, double[] latlng, double[] lunarAscDecDist) {
        double jce = calcTimeJulianCent(jd);
        double sidereal = greenwichMeanSiderealTime(jce);
        double lunarDistance = lunarAscDecDist[coords.DISTANCE.ordinal()];
        double standardAltitude = degToRad(-.583 - equatorialHorizontalParallax(lunarDistance));

        double lat = latlng[coords.LAT.ordinal()];
        double lng = -1 * latlng[coords.LNG.ordinal()]; // reverse the longitude, so it is measured positively west from Greenwich.

        double localHourAngle = Math.sin(standardAltitude) - (Math.sin(degToRad(lat)) * Math.sin(degToRad(lunarAscDecDist[position.DECLINATION.ordinal()]))) / (Math.cos(degToRad(lat)) * Math.cos(degToRad(lunarAscDecDist[position.DECLINATION.ordinal()])));
        localHourAngle = radToDeg(Math.acos(localHourAngle));

        double transit = getInRange((lunarAscDecDist[position.ASCENSION.ordinal()] + lng - sidereal) / 360);
        double rising = transit - (localHourAngle / 360);
        double setting = transit + (localHourAngle / 360);

        return new double[]{rising, transit, setting};
    }

    /**
     * Calculates precise times for the Moon's rising, transit, and setting at a given location
     * and date using an iterative refinement algorithm. This method improves upon the initial
     * approximations by repeatedly recalculating the Moon's position at the estimated event times.

     * The Moon moves approximately 13° per day relative to the stars, which is about 0.5° per hour.
     * This rapid motion means that approximate calculations can be off by several minutes. This
     * method performs six iterations to converge on accurate event times, accounting for the
     * Moon's changing position during the day.

     * The algorithm handles edge cases where events occur on different days (e.g., when the
     * approximate time falls before midnight or after the next midnight), ensuring the returned
     * times are for events near the specified date.
     *
     * @param lat Observer's latitude in degrees (positive north, negative south)
     * @param lng Observer's longitude in degrees (positive east, negative west)
     * @param date The date for which to calculate rise, transit, and set times
     * @return List of ZonedDateTime objects in order: [moonrise, transit, moonset], converted
     *         to the system's local time zone. Returns times in UTC if events don't occur on
     *         the specified date (circumpolar or never rises).
     */
    public List<ZonedDateTime> moonRisingSettingTransitPrecise(double lat, double lng, LocalDate date) {
        List<ZonedDateTime> zonedDateTimeList = new ArrayList<>();
        int[] positions = {time.Rise.ordinal(), time.Transit.ordinal(), time.Set.ordinal()};
        double[] coords = {lat, lng};

        for (int position : positions
        ) {
            double timeDifference = 0;
            LocalDate tempDate = date;
            for (int i = 0; i < 6; i++) {
                double jd = getJDFromCalenderDate(tempDate.getYear(), tempDate.getMonthValue(), tempDate.getDayOfMonth() + timeDifference);
                double[] td = lunarAscensionDeclinationDistance(jd);
                double[] riseTransitSet = risingTransitSettingApproximate(jd, coords, td);
                double fractionOfDay = riseTransitSet[position];
                if (fractionOfDay < 0) {
                    tempDate.minusDays((long) Math.abs(Math.floor(fractionOfDay)));
                    timeDifference = 1.0 + (Math.abs(Math.ceil(fractionOfDay)) + fractionOfDay);
                } else if (fractionOfDay > 1) {
                    tempDate.plusDays((long) Math.floor(fractionOfDay));
                    timeDifference = (fractionOfDay - Math.floor(fractionOfDay));
                } else timeDifference = fractionOfDay;
            }
            zonedDateTimeList.add(position, localFractionOfDayFromUTCToLocal(timeDifference, tempDate));
        }
        return zonedDateTimeList;
    }
}