plugins {
    id("java")
    id("maven-publish")
}

group = "io.github.dgrims3"
version = "1.0.0"

repositories {
    mavenCentral()
}

java {
    toolchain {
        // 1. Run the build using JDK 21 (your machine's version)
        languageVersion.set(JavaLanguageVersion.of(21))
    }
    // Generate sources and javadoc JARs for publication
    withSourcesJar()
    withJavadocJar()
}

tasks.withType<JavaCompile>().configureEach {
    // 2. Force the output to be compatible with Java 17
    options.release.set(17)
}

tasks.javadoc {
    if (JavaVersion.current().isJava9Compatible) {
        (options as StandardJavadocDocletOptions).addBooleanOption("html5", true)
    }
}

dependencies {
    testImplementation(platform("org.junit:junit-bom:5.9.1"))
    testImplementation("org.junit.jupiter:junit-jupiter")
}

tasks.test {
    useJUnitPlatform()
}

publishing {
    repositories {
        maven {
            name = "GitHubPackages"
            url = uri("https://maven.pkg.github.com/dgrims3/MoonLocatorLibrary")
            credentials {
                username = project.findProperty("gpr.user") as String? ?: System.getenv("GITHUB_ACTOR")
                password = project.findProperty("gpr.token") as String? ?: System.getenv("GITHUB_TOKEN")
            }
        }
    }
    publications {
        create<MavenPublication>("gpr") {
            from(components["java"])

            pom {
                name.set("MoonLocatorLibrary")
                description.set("A Java library for calculating moon positions, rise/set times, and lunar properties")
                url.set("https://github.com/dgrims3/MoonLocatorLibrary")

                licenses {
                    license {
                        name.set("MIT License")
                        url.set("https://opensource.org/licenses/MIT")
                    }
                }

                developers {
                    developer {
                        id.set("dgrims3")
                        name.set("David Grimes")
                        email.set("dgrims3@github.com")
                    }
                }

                scm {
                    connection.set("scm:git:git://github.com/dgrims3/MoonLocatorLibrary.git")
                    developerConnection.set("scm:git:ssh://github.com/dgrims3/MoonLocatorLibrary.git")
                    url.set("https://github.com/dgrims3/MoonLocatorLibrary")
                }
            }
        }
    }
}
