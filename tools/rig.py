#!/usr/bin/env python3
"""
Shared bits every sweep needs: where the compiled classes and the dependency jars are.

The classpath is not checked in - it is a list of absolute paths into ~/.m2 that differ between
machines. Generate it once with

    mvn -o dependency:build-classpath -Dmdep.outputFile=tools/cp.txt

or point SMFRET_CP at one.
"""
import os
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

CLASSES = os.path.join(ROOT, "target", "classes")


def _jdk(tool):
    """SMFRET_JAVA(C), then JAVA_HOME, then PATH.

    Unlike the Maven build these do not need JDK 11 - that requirement comes from pom-scijava
    asking for Java 8 bytecode, and nothing here is compiled with --release. Any modern JDK will
    do, so falling back to whatever is on PATH is safe.
    """
    override = os.environ.get("SMFRET_" + tool.upper())
    if override:
        return override

    home = os.environ.get("JAVA_HOME")
    if home:
        candidate = os.path.join(home, "bin", tool)
        if os.path.exists(candidate):
            return candidate

    found = shutil.which(tool)
    if found:
        return found
    raise SystemExit(f"No {tool}. Set JAVA_HOME or SMFRET_{tool.upper()}.")


JAVA = _jdk("java")
JAVAC = _jdk("javac")


def classpath():
    """The dependency classpath, from SMFRET_CP or tools/cp.txt."""
    direct = os.environ.get("SMFRET_CP")
    if direct:
        return direct

    path = os.path.join(HERE, "cp.txt")
    if not os.path.exists(path):
        raise SystemExit(
            "No classpath. Run:\n"
            "  mvn -o dependency:build-classpath -Dmdep.outputFile=tools/cp.txt\n"
            "or set SMFRET_CP.")
    return open(path).read().strip()


def compile_harnesses(*names):
    """Compile the Java drivers in this directory into the working directory."""
    import subprocess
    sources = [os.path.join(HERE, n) for n in names]
    subprocess.run([JAVAC, "-nowarn", "-cp", os.pathsep.join([CLASSES, classpath()]),
                    "-d", ".", *sources], check=True)
