"""Una corrida: tomar, preparar, lanzar, vigilar, post-procesar, registrar.

La falla de una corrida nunca corta el programa: queda escrita en su archivo `status`.
"""
import glob
import os
import re
import shutil
import signal
import socket
import subprocess
import time
import traceback

import cfg
import parsers
import post

TERMINAL = ("ok", "sin_converger", "divergio", "fallo", "vencida", "estancada", "interrumpida")


def now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


# ---------------------------------------------------------------- estado

def write_status(rd, state, reason=""):
    reason = " ".join(str(reason).split())[:300]
    tmp = os.path.join(rd, "status.tmp")
    with open(tmp, "w", encoding="utf8") as f:
        f.write("%s | %s | %s | %s\n" % (state, now(), socket.gethostname(), reason))
    os.replace(tmp, os.path.join(rd, "status"))


def read_status(rd):
    try:
        with open(os.path.join(rd, "status"), encoding="utf8") as f:
            parts = [p.strip() for p in f.readline().split("|", 3)]
    except OSError:
        parts = ["pendiente"]
    parts += [""] * (4 - len(parts))
    return dict(zip(("estado", "fecha", "host", "motivo"), parts))


def claim(rd):
    try:
        fd = os.open(os.path.join(rd, "lock"), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as f:
        f.write("%s %d %s\n" % (socket.gethostname(), os.getpid(), now()))
    return True


def release(rd):
    try:
        os.remove(os.path.join(rd, "lock"))
    except OSError:
        pass


def heartbeat_age(rd):
    try:
        return time.time() - os.path.getmtime(os.path.join(rd, "lock"))
    except OSError:
        return None


# ---------------------------------------------------------------- vigilancia

class Scanner:
    """Lee sólo lo agregado desde la última llamada y busca patrones de falla."""

    def __init__(self, patterns):
        self.patterns = [(state, re.compile(rx, re.M)) for state, rx in patterns]
        self.offset = {}

    def scan(self, paths, final=False):
        for path in paths:
            try:
                size = os.path.getsize(path)
            except OSError:
                continue
            start = self.offset.get(path, 0)
            if size < start:
                start = 0
            with open(path, "rb") as f:
                f.seek(start)
                chunk = f.read(size - start)
            end = len(chunk) if final else chunk.rfind(b"\n") + 1
            if end <= 0:
                continue
            self.offset[path] = start + end
            text = chunk[:end].decode("utf8", "replace")
            for state, rx in self.patterns:
                m = rx.search(text)
                if m:
                    a = text.rfind("\n", 0, m.start()) + 1
                    b = text.find("\n", m.end())
                    return state, "%s: %s" % (os.path.basename(path), text[a:b if b >= 0 else None].strip())
        return None


def patterns(caso):
    extra = [(state, caso.get("watch", state)) for state in ("divergio", "fallo")]
    return parsers.PATTERNS.get(caso.solver, []) + [(s, rx) for s, rx in extra if rx]


def expand(rd, globs):
    out = []
    for g in globs:
        out.extend(sorted(glob.glob(os.path.join(rd, g))))
    return out


def kill(p):
    for sig, wait in ((signal.SIGTERM, 20), (signal.SIGKILL, 5)):
        try:
            os.killpg(p.pid, sig)
        except (ProcessLookupError, PermissionError):
            return
        try:
            p.wait(wait)
            return
        except subprocess.TimeoutExpired:
            pass


def shell(script, mach, solver):
    env = "\n".join(mach.get(k, "") for k in ("env", "env_" + solver))
    return ["bash", "-c", (env + "\n" + script).strip()]


def supervise(argv, cwd, logpath, rd, limits, watch=(), scanner=None):
    """Corre y vigila. Devuelve (código de salida, veredicto o None)."""
    timeout, stall, poll = limits
    with open(logpath, "ab") as lf:
        p = subprocess.Popen(argv, cwd=cwd, stdout=lf, stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL, start_new_session=True)
    t0 = last_change = time.time()
    size, verdict = -1, None
    try:
        while verdict is None:
            try:
                p.wait(poll)
                break
            except subprocess.TimeoutExpired:
                pass
            try:
                os.utime(os.path.join(rd, "lock"))
            except OSError:
                pass
            t = time.time()
            paths = [logpath] + expand(rd, watch)
            current = sum(os.path.getsize(x) for x in paths if os.path.exists(x))
            if current != size:
                size, last_change = current, t
            if timeout and t - t0 > timeout:
                verdict = ("vencida", "superó el límite de %g s" % timeout)
            elif stall and t - last_change > stall:
                verdict = ("estancada", "sin salida nueva durante %g s" % stall)
            elif scanner:
                verdict = scanner.scan(paths)
    except KeyboardInterrupt:
        kill(p)
        raise
    if verdict:
        kill(p)
    return p.wait(), verdict


def tail(path, n=1):
    try:
        with open(path, "rb") as f:
            f.seek(max(0, os.path.getsize(path) - 4096))
            lines = f.read().decode("utf8", "replace").strip().splitlines()
        return " / ".join(lines[-n:])
    except OSError:
        return ""


# ---------------------------------------------------------------- corrida

def execute(caso, row, mach, root=None):
    rd = caso.run_dir(row["id"], mach, root)
    os.makedirs(rd, exist_ok=True)
    if not claim(rd):
        return "ocupada"
    try:
        state, _ = _execute(caso, row, mach, rd)
    except KeyboardInterrupt:
        write_status(rd, "interrumpida", "señal de interrupción")
        raise
    except Exception:
        with open(os.path.join(rd, "runner.err"), "w", encoding="utf8") as f:
            f.write(traceback.format_exc())
        state = "fallo"
        write_status(rd, state, "runner: %s (runner.err)" % traceback.format_exc().strip().splitlines()[-1])
    finally:
        release(rd)
    return state


def limits(caso):
    get = lambda k, d: float(caso.get("run", k, d))
    return get("timeout_s", "0"), get("stall_s", "0"), get("poll_s", "30")


def _execute(caso, row, mach, rd):
    case = os.path.join(rd, "case")
    for name in ("log", "result.ini", "prepare.log", "runner.err"):
        if os.path.exists(os.path.join(rd, name)):
            os.remove(os.path.join(rd, name))
    if os.path.exists(case):
        shutil.rmtree(case + ".anterior", ignore_errors=True)
        os.replace(case, case + ".anterior")

    try:
        vals = cfg.resolve(caso, row, mach)
        prepare = cfg.template(caso, "prepare", vals, rd)
        command = cfg.template(caso, "command", vals, rd)
        watch = (cfg.template(caso, "watch", vals, rd) or "").split()
        if not command:
            raise ValueError("[run] command vacío")
    except ValueError as e:
        write_status(rd, "fallo", "configuración: %s" % e)
        return "fallo", str(e)

    meta = {"id": row["id"], "caso": caso.name, "maquina": mach["nombre"],
            "host": socket.gethostname(), "fecha": now(), "taubench_s": mach["taubench_s"]}
    cfg.write_ini(os.path.join(rd, "run.ini"),
                  {"corrida": meta, "valores": vals,
                   "comandos": {"prepare": prepare or "", "command": command}})
    timeout, stall, poll = limits(caso)

    if prepare:
        write_status(rd, "corriendo", "preparando")
        plog = os.path.join(rd, "prepare.log")
        rc, verdict = supervise(shell(prepare, mach, caso.solver), rd, plog, rd, (3600, 0, min(poll, 5)))
        if rc != 0 or verdict:
            why = verdict[1] if verdict else "código %d" % rc
            write_status(rd, "fallo", "preparación: %s; %s" % (why, tail(plog)))
            return "fallo", why
    if not os.path.isdir(case):
        if prepare:
            write_status(rd, "fallo", "preparación no creó case/")
            return "fallo", "sin case/"
        os.makedirs(case)

    write_status(rd, "corriendo", "")
    log = os.path.join(rd, "log")
    scanner = Scanner(patterns(caso))
    rc, verdict = supervise(shell(command, mach, caso.solver), case, log, rd, (timeout, stall, poll), watch, scanner)
    verdict = verdict or scanner.scan([log] + expand(rd, watch), final=True)

    write_status(rd, "corriendo", "post-proceso")
    return finish(caso, rd, vals, meta, rc, verdict)


def finish(caso, rd, vals, meta, rc=0, verdict=None):
    res, perr = None, ""
    try:
        res = post.process(caso, rd, vals, meta)
        post.write(rd, res)
    except Exception as e:
        res, perr = None, "%s: %s" % (type(e).__name__, e)
        with open(os.path.join(rd, "runner.err"), "w", encoding="utf8") as f:
            f.write(traceback.format_exc())

    if verdict:
        state, reason = verdict
    elif rc != 0:
        state, reason = "fallo", "código de salida %d; %s" % (rc, tail(os.path.join(rd, "log")))
    elif res is None:
        state, reason = "fallo", "post-proceso: " + perr
    elif post.kind(caso) == "errores":
        state, reason = ("ok", "") if res["errores"] else ("fallo", "el log no tiene normas de error")
    elif not res["historia"]["filas"]:
        state, reason = "fallo", "terminó sin historial de fuerzas"
    elif res["parada"]["cumplio"] == "si":
        state, reason = "ok", ""
    else:
        state, reason = "sin_converger", res["parada"].get("nota", "")
    write_status(rd, state, reason)
    return state, reason


def repost(caso, rd):
    """Rehace result.ini con la configuración actual de [post]. Sólo reclasifica ok/sin_converger."""
    run = cfg.read_ini(os.path.join(rd, "run.ini"))
    vals = {k: cfg.number(v) for k, v in run["valores"].items()}
    source = "log" if post.kind(caso) == "errores" else "forces"
    if not glob.glob(cfg.template(caso, source, vals, rd, "post") or ""):
        return read_status(rd)["estado"]
    meta = dict(run["corrida"])
    meta["taubench_s"] = cfg.taubench_of(meta.get("maquina")) or cfg.number(meta.get("taubench_s") or "") or None
    before = read_status(rd)
    if before["estado"] in ("ok", "sin_converger"):
        return finish(caso, rd, vals, meta)[0]
    try:
        post.write(rd, post.process(caso, rd, vals, meta))
    except Exception:
        pass
    return before["estado"]
