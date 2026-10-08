#!/usr/bin/env python3
"""PreToolUse: protege base (/opt/nektar++, ~/nektar/build), raw/ y master. Silencioso si permite."""
import json, os, re, subprocess, sys

H = os.path.expanduser("~")
NEK = f"{H}/nektar"
TESIS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def deny(msg):
    print(f"guard: {msg}", file=sys.stderr)
    sys.exit(2)

def branch():
    try:
        return subprocess.run(["git", "-C", NEK, "branch", "--show-current"],
                              capture_output=True, text=True, timeout=5).stdout.strip()
    except Exception:
        return ""

ev = json.load(sys.stdin)
tool, inp = ev.get("tool_name", ""), ev.get("tool_input", {})

if tool in ("Edit", "Write", "NotebookEdit"):
    p = os.path.realpath(os.path.expanduser(inp.get("file_path") or inp.get("notebook_path") or ""))
    if p.startswith("/opt/nektar++"):
        deny("/opt/nektar++ es la instalación base; no se edita.")
    if p.startswith(f"{NEK}/build/") or p.startswith(f"{NEK}/build-svv/"):
        deny("los directorios de build no se editan a mano.")
    if p.startswith(f"{TESIS}/raw/"):
        deny("raw/ contiene fuentes; sólo lectura.")
    if p.startswith(NEK + "/"):
        b = branch()
        if b != "feature/svv-cfs":
            deny(f"~/nektar está en '{b}'; editar sólo en feature/svv-cfs (git -C ~/nektar switch feature/svv-cfs).")
        if not p.startswith(f"{NEK}/solvers/CompressibleFlowSolver/"):
            deny("fuera de solvers/CompressibleFlowSolver/; registrar la excepción en DECISIONES.md y pedir confirmación al usuario.")

elif tool == "Bash":
    c = inp.get("command", "")
    rules = [
        (r"\bsudo\b", "sin sudo."),
        (r"git\b[^|;&]*\bpush\b", "git push no permitido desde agentes."),
        (r"(rm|mv)\s+(-\w+\s+)*[^|;&]*(/opt/nektar\+\+|~/nektar/?(\s|$)|" + re.escape(NEK) + r"/?(\s|$)|nektar/build(\s|/|$))",
         "no se borra ni mueve la instalación base, el repo ni su build."),
        (r"\b(make|ninja)\b[^|;&]*\binstall\b", "instalar sólo vía harness/bin/svv-build (prefijo build-svv/dist)."),
        (r"cmake\s+--install(?![^|;&]*build-svv)", "cmake --install sólo sobre build-svv."),
        (r"cmake\s+--build\s+\S*nektar/build(\s|$)", "~/nektar/build es la base; usar harness/bin/svv-build."),
        (r"git\s+(-C\s+\S*nektar\s+)?(reset\s+--hard|clean\s+-\w*f|checkout\s+--\s|restore\s)", "operación git destructiva; pedir confirmación al usuario."),
    ]
    for pat, msg in rules:
        if re.search(pat, c):
            deny(msg)
sys.exit(0)
