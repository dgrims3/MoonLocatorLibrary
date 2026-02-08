# MoonLocator

A high-precision Java library for calculating lunar positions, phases, rise/set times, and other moon-related astronomical data. Built on astronomical algorithms from NREL and Jean Meeus.

## Origin

This library was originally created for the Android application [Lunar Locator](https://play.google.com/store/apps/details?id=com.ll.LunarLocator), available on Google Play Store. It has since been extracted and published as a standalone library for broader use in astronomical calculations and celestial navigation applications.

## Features

- **Lunar Position Calculations**
  - Azimuth and altitude (horizontal coordinates)
  - Right ascension and declination (equatorial coordinates)
  - Geographic latitude and longitude of moon's sub-point
  - Local hour angle

- **Rise/Set/Transit Times**
  - Precise moonrise times
  - Moonset times
  - Moon transit (culmination) times
  - Iterative refinement for high accuracy

- **Lunar Properties**
  - Earth-Moon distance in kilometers
  - Illumination percentage (lunar phase)
  - Equatorial horizontal parallax

- **Flexible Time Zone Support**
  - Works with `ZonedDateTime` and `LocalDateTime`
  - Automatic UTC offset handling
  - Configurable time zones

## Installation

### Maven Central (Recommended)

Add the following to your `build.gradle.kts`:

```kotlin
repositories {
    mavenCentral()
}

dependencies {
    implementation("io.github.dgrims3:moon-library:1.0.0")
}
```

Or for Maven, add to your `pom.xml`:

```xml
<dependency>
    <groupId>io.github.dgrims3</groupId>
    <artifactId>moon-library</artifactId>
    <version>1.0.0</version>
</dependency>
```

### JitPack

Alternatively, you can use JitPack:

```kotlin
repositories {
    mavenCentral()
    maven { url = uri("https://jitpack.io") }
}

dependencies {
    implementation("com.github.dgrims3:MoonLocatorLibrary:v1.0.0")
}
```

### Requirements

- Java 17 or higher (compiled with Java 21, compatible with Java 17+)

## Usage

### Basic Example

```java
import io.github.dgrims3.MoonLocator;
import java.time.LocalDateTime;
import java.time.ZonedDateTime;
import java.time.ZoneId;

// Create a MoonLocator for a specific location and time
LocalDateTime dateTime = LocalDateTime.of(2025, 1, 15, 20, 30);
double latitude = 40.7128;  // New York City latitude
double longitude = -74.0060; // New York City longitude

MoonLocator locator = new MoonLocator(dateTime, latitude, longitude);

// Get moon position
double azimuth = locator.getAzimuth();
double altitude = locator.getAltitude();
System.out.println("Moon azimuth: " + azimuth + "°");
System.out.println("Moon altitude: " + altitude + "°");

// Get rise/set times
ZonedDateTime moonrise = locator.getMoonRise();
ZonedDateTime moonset = locator.getMoonSet();
ZonedDateTime transit = locator.getMoonTransit();
System.out.println("Moonrise: " + moonrise);
System.out.println("Moonset: " + moonset);
System.out.println("Transit: " + transit);

// Get lunar properties
double distance = locator.getDistance();
double illumination = locator.getIllumination();
System.out.println("Distance: " + distance + " km");
System.out.println("Illumination: " + illumination + "%");
```

### With Specific Time Zone

```java
ZoneId zoneId = ZoneId.of("America/New_York");
MoonLocator locator = new MoonLocator(dateTime, latitude, longitude, zoneId);
```

### With ZonedDateTime

```java
ZonedDateTime zonedDateTime = ZonedDateTime.now(ZoneId.of("UTC"));
MoonLocator locator = new MoonLocator(zonedDateTime, latitude, longitude);
```

### Geographic Sub-Point

```java
// Get the point on Earth where the moon is directly overhead
double moonLat = locator.getMoonGeographicLatitude();
double moonLng = locator.getMoonGeographicLongitude();
System.out.println("Moon is directly overhead at: " + moonLat + ", " + moonLng);
```

### Equatorial Coordinates

```java
double rightAscension = locator.getRightAscension();
double declination = locator.getDeclination();
double hourAngle = locator.getLocalHourAngle();
```

## API Reference

### Constructor

- `MoonLocator(ZonedDateTime dateTime, double latitude, double longitude)`
- `MoonLocator(LocalDateTime dateTime, double latitude, double longitude)`
- `MoonLocator(LocalDateTime dateTime, double latitude, double longitude, ZoneId zoneId)`

**Parameters:**
- `dateTime` - The date and time for calculations
- `latitude` - Observer's latitude in degrees (positive north, negative south)
- `longitude` - Observer's longitude in degrees (positive east, negative west)
- `zoneId` - Time zone for the calculation (optional, defaults to system time zone)

### Position Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `getAzimuth()` | `double` | Compass direction from north (0-360°) |
| `getAltitude()` | `double` | Angle above horizon (-90 to +90°) |
| `getRightAscension()` | `double` | Right ascension in degrees |
| `getDeclination()` | `double` | Declination in degrees |
| `getLocalHourAngle()` | `double` | Local hour angle in degrees |

### Time Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `getMoonRise()` | `ZonedDateTime` | Time when moon crosses horizon ascending |
| `getMoonSet()` | `ZonedDateTime` | Time when moon crosses horizon descending |
| `getMoonTransit()` | `ZonedDateTime` | Time when moon reaches highest point |

### Property Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `getDistance()` | `double` | Earth-Moon distance in kilometers |
| `getIllumination()` | `double` | Illuminated percentage (0-100) |
| `getMoonGeographicLatitude()` | `double` | Latitude of moon's sub-point |
| `getMoonGeographicLongitude()` | `double` | Longitude of moon's sub-point |

### Accessor Methods

| Method | Return Type | Description |
|--------|-------------|-------------|
| `getDateTime()` | `ZonedDateTime` | Configured date and time |
| `getZoneId()` | `ZoneId` | Configured time zone |
| `getLatitude()` | `double` | Observer's latitude |
| `getLongitude()` | `double` | Observer's longitude |
| `getOffSet()` | `int` | UTC offset in hours |

## Coordinate Systems

### Horizontal Coordinates (Alt-Az)
- **Azimuth**: 0° = North, 90° = East, 180° = South, 270° = West
- **Altitude**: 0° = Horizon, 90° = Zenith, negative = below horizon

### Equatorial Coordinates
- **Right Ascension**: Celestial longitude (0-360°)
- **Declination**: Celestial latitude (+90° to -90°)

### Geographic Coordinates
- **Latitude**: Positive north, negative south
- **Longitude**: Positive east, negative west

## Algorithm References

This library implements algorithms from:

1. **Reda, I. (2010)**. Solar Eclipse Monitoring for Solar Energy Applications Using the Solar and Moon Position Algorithms. NREL/TP-3B0-47681. National Renewable Energy Laboratory.

2. **Meeus, J. (1998)**. Astronomical Algorithms. 2nd Edition. Willmann-Bell, Inc.

The implementation uses:
- Chapront's ELP-2000/82 lunar theory for position calculations
- Periodic term summations with 60 terms for longitude/distance and latitude
- Iterative refinement for rise/set/transit times
- Atmospheric refraction and parallax corrections

## Accuracy

- Position accuracy: Sub-arcminute for dates within ±100 years of J2000.0
- Rise/set times: Within ~1 minute for most locations
- Distance: Within ~1 kilometer
- Illumination: Within 0.1%

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Author

**David Grimsley** ([@dgrims3](https://github.com/dgrims3))

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Related Projects

- [Lunar Locator](https://play.google.com/store/apps/details?id=com.ll.LunarLocator) - Android application using this library

## Support

For issues, questions, or contributions, please visit the [GitHub repository](https://github.com/dgrims3/MoonLocatorLibrary).
