"""Pruebas: parsers contra salidas reales, numérica, runner ante cada falla, reportes.

    cd casos && python3 -m unittest discover -s tests -v
"""
import csv
import datetime
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import textwrap
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import cfg        # noqa: E402
import num        # noqa: E402
import parsers    # noqa: E402
import report     # noqa: E402
import runner     # noqa: E402
import telemetry  # noqa: E402


def fix(name):
    return os.path.join(HERE, "fixtures", name)


def data_rows(path):
    with open(path) as f:
        return sum(1 for line in f if line.strip() and not line.lstrip().startswith("#"))


def scan(solver, path):
    return runner.Scanner(parsers.PATTERNS[solver]).scan([path], final=True)


def caso_de_prueba(tmp, text):
    d = os.path.join(tmp, "prueba")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "caso.ini"), "w") as f:
        f.write(textwrap.dedent(text))
    caso = cfg.Caso(d)
    if caso.ini.has_section("sweep"):
        cfg.write_matrix(caso.matrix_path, cfg.plan(caso))
    return caso


def machine(tmp, taubench=None):
    return {"nombre": "prueba", "procs": 2, "runs_root": os.path.join(tmp, "runs"),
            "taubench_s": taubench, "vars": {}}


# ---------------------------------------------------------------- parsers

class Parsers(unittest.TestCase):
    def test_nektar_forces(self):
        f = parsers.forces("nektar", [fix("nektar.fce")])
        self.assertEqual(len(f["t"]), data_rows(fix("nektar.fce")))
        self.assertEqual(f["t"][0], 0.0)
        self.assertEqual(len(f["t"]), len(f["fy"]))

    def test_nektar_timing_discounts_writes(self):
        t, s = parsers.timing("nektar", fix("nektar.log"))
        self.assertEqual(len(t), 5)
        self.assertTrue(all(b >= a for a, b in zip(s, s[1:])))
        with open(fix("nektar.log")) as f:
            total = float(re.search(r"Time-integration\s*:\s*(\S+)s", f.read()).group(1))
        self.assertLess(s[-1], total)
        self.assertGreater(s[-1], 0.5 * total)

    def test_openfoam(self):
        path = fix("openfoam_force.dat")
        self.assertEqual(len(parsers.forces("openfoam", [path])["t"]), data_rows(path))
        t, s = parsers.timing("openfoam", fix("openfoam.log"))
        self.assertTrue(t)
        self.assertEqual(len(t), len(s))

    def test_failure_signatures(self):
        for solver, name in (("nektar", "nektar.log"), ("nektar", "nektar.fce"),
                             ("openfoam", "openfoam.log"), ("openfoam", "openfoam_force.dat")):
            self.assertIsNone(scan(solver, fix(name)), name)
        self.assertEqual(scan("nektar", fix("nektar_nan.log"))[0], "divergio")
        self.assertEqual(scan("openfoam", fix("openfoam_fail.log"))[0], "divergio")

    def test_nan_row_in_forces(self):
        with tempfile.NamedTemporaryFile("w", suffix=".fce", delete=False) as f:
            f.write("#   Time  F1-total\n 0.1  1.0\n 0.2  -nan\n")
        self.assertEqual(scan("nektar", f.name)[0], "divergio")
        os.remove(f.name)

    def test_error_norms(self):
        e = parsers.errors(fix("nektar.log"))
        with open(fix("nektar.log")) as f:
            text = f.read()
        expected = float(re.search(r"L inf error \(variable E\) : (\S+)", text).group(1))
        self.assertEqual(e["Linf_E"], expected)
        self.assertEqual(len(e), 8)

    def test_count_cells(self):
        mesh = os.path.join(ROOT, "..", "nektar_cases", "nektar_paper", "naca0012", "naca0012.xml")
        if not os.path.exists(mesh):
            self.skipTest("sin malla NACA")
        self.assertEqual(parsers.count_cells("nektar", mesh), 3002)


# ---------------------------------------------------------------- numérica

class Numerics(unittest.TestCase):
    def test_richardson(self):
        q = lambda h: 1.0 + 0.5 * h**2
        rich, _ = num.richardson(q(1), q(0.5), q(0.25))
        self.assertAlmostEqual(rich["p"], 2.0, 9)
        self.assertAlmostEqual(rich["qinf"], 1.0, 9)
        self.assertIsNone(num.richardson(1.0, 0.9, 0.95)[0])

    def test_pchip(self):
        xs, ys = [0, 1, 2, 3], [0, 0, 1, 1]
        values = [num.pchip(xs, ys, 3 * i / 100) for i in range(101)]
        self.assertTrue(all(-1e-12 <= v <= 1 + 1e-12 for v in values))
        self.assertAlmostEqual(num.pchip([0, 1, 2], [0, 2, 4], 1.5), 3.0)
        self.assertIsNone(num.pchip(xs, ys, 3.5))

    def test_first_settled(self):
        x = [0.1 * i for i in range(200)]
        s = [math.exp(-v) for v in x]
        j, lo = num.first_settled(x, {"cd": s}, 1.0, {"cd": 1e-3})
        self.assertLess(s[lo] - s[j], 1e-3)
        lo2 = num.window_start(x, j - 1, 1.0)
        self.assertGreaterEqual(s[lo2] - s[j - 1], 1e-3)


# ---------------------------------------------------------------- runner

EXPECTED = {"ok": "ok", "noconv": "sin_converger", "nan": "divergio", "crash": "fallo",
            "fatal0": "fallo", "hang": "estancada", "slow": "vencida", "badcfg": "fallo"}

RUNNER_INI = """
[caso]
solver = nektar
id     = {modo}

[sweep]
modo = ok noconv nan crash fatal0 hang slow badcfg

[values]
rho  = 2.0
U    = 1.0 if modo != 'badcfg' else no_existe
aref = 1.0
lref = 1.0

[run]
command   = python3 {tools}/tests/fake_solver.py {modo}
watch     = case/forces.fce
timeout_s = 4
stall_s   = 2
poll_s    = 0.2

[post]
forces = {case}/forces.fce
log    = {run}/log
unit   = conv
window = 0.05
tol_cd = 1e-4

[report]
kind = cost_curve
"""


class Runner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="caso_")
        cls.caso = caso_de_prueba(cls.tmp, RUNNER_INI)
        cls.mach = machine(cls.tmp, taubench=4.0)
        cls.states = {r["id"]: runner.execute(cls.caso, r, cls.mach) for r in cls.caso.rows()}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def rd(self, rid):
        return self.caso.run_dir(rid, self.mach)

    def test_each_failure_is_classified(self):
        self.assertEqual(self.states, EXPECTED)
        for rid, state in EXPECTED.items():
            self.assertEqual(runner.read_status(self.rd(rid))["estado"], state)
            self.assertFalse(os.path.exists(os.path.join(self.rd(rid), "lock")))

    def test_ok_result_is_normalized(self):
        r = cfg.read_ini_flat(os.path.join(self.rd("ok"), "result.ini"))
        self.assertEqual(r["parada.cumplio"], "si")
        self.assertAlmostEqual(float(r["resultado.cd"]), 0.02, delta=1e-4)
        self.assertAlmostEqual(float(r["resultado.cl"]), 0.3, 9)
        s = float(r["costo.segundos_solver"])
        self.assertGreater(s, 0)
        self.assertAlmostEqual(float(r["costo.costo_normalizado"]), 2 * s / 4.0, 6)

    def test_failed_runs_keep_partial_result(self):
        r = cfg.read_ini_flat(os.path.join(self.rd("nan"), "result.ini"))
        self.assertGreater(int(r["historia.filas"]), 0)

    def test_no_orphan_processes(self):
        out = subprocess.run(["pgrep", "-f", "fake_solver.py"], capture_output=True, text=True).stdout
        self.assertEqual(out.strip(), "")

    def test_report_and_rerun(self):
        text = report.regenerate(self.caso, [self.mach["runs_root"]])
        self.assertIn("ok 1", text)
        self.assertIn("divergio 1", text)
        with open(os.path.join(self.caso.results_dir, "results.csv")) as f:
            self.assertEqual(len(f.readlines()), 1 + len(EXPECTED))
        row = [r for r in self.caso.rows() if r["id"] == "ok"][0]
        self.assertEqual(runner.execute(self.caso, row, self.mach), "ok")
        self.assertTrue(os.path.isdir(os.path.join(self.rd("ok"), "case.anterior")))

    def test_repost_changes_classification(self):
        path = os.path.join(self.caso.dir, "caso.ini")
        with open(path) as f:
            original = f.read()
        try:
            with open(path, "w") as f:
                f.write(original.replace("tol_cd = 1e-4", "tol_cd = 1e-12"))
            caso = cfg.Caso(self.caso.dir)
            self.assertEqual(runner.repost(caso, self.rd("ok")), "sin_converger")
            self.assertEqual(runner.repost(caso, self.rd("nan")), "divergio")
        finally:
            with open(path, "w") as f:
                f.write(original)
            runner.repost(cfg.Caso(self.caso.dir), self.rd("ok"))


class Interrupt(unittest.TestCase):
    def test_sigint_marks_and_cleans(self):
        tmp = tempfile.mkdtemp(prefix="caso_")
        caso = caso_de_prueba(tmp, RUNNER_INI.replace("ok noconv nan crash fatal0 hang slow badcfg", "hang")
                        .replace("timeout_s = 4", "timeout_s = 0").replace("stall_s   = 2", "stall_s = 0"))
        rd = caso.run_dir("hang", machine(tmp))
        script = ("import sys; sys.path.insert(0, %r)\nimport cfg, runner\n"
                  "caso = cfg.Caso(%r)\n"
                  "try:\n    runner.execute(caso, caso.rows()[0], %r)\n"
                  "except KeyboardInterrupt:\n    pass\n") % (ROOT, caso.dir, machine(tmp))
        p = subprocess.Popen([sys.executable, "-c", script])
        for _ in range(100):
            if os.path.exists(os.path.join(rd, "case", "forces.fce")):
                break
            time.sleep(0.1)
        time.sleep(0.5)
        p.send_signal(signal.SIGINT)
        p.wait(40)
        self.assertEqual(runner.read_status(rd)["estado"], "interrumpida")
        self.assertFalse(os.path.exists(os.path.join(rd, "lock")))
        out = subprocess.run(["pgrep", "-f", "fake_solver.py hang"], capture_output=True, text=True).stdout
        self.assertEqual(out.strip(), "")
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------- reportes y telemetría

def synthetic_flight(path, cd=0.45, m=20.0, A=0.0179):
    """Vuelo balístico con Cd constante, 10 Hz, en el formato de 16 columnas del archivo real."""
    v, h, g, t, dt = 420.0, 1500.0, math.radians(80), 0.0, 1e-3
    base = datetime.datetime(2026, 6, 17, 11, 0, 0)
    with open(path, "w") as f:
        f.write("hora,idx,tipo,lat,lon,alt,a,b,c,d,e,vel,f,g,h,i\n")
        f.write("%s,0,OTHER,0,0,0,0,0,0,0,0,0,0,0,0,0\n" % base.isoformat())
        for k in range(120001):
            rho, _ = telemetry.isa(h)
            drag = 0.5 * rho * v * v * cd * A
            if k % 100 == 0:
                stamp = (base + datetime.timedelta(seconds=t)).isoformat()
                f.write("%s,%d,TRACK,31,-103,%.6f,3,%.8f,0,0,0,%.6f,0,0,0,1\n"
                        % (stamp, k, h, -drag / m, v))
            v, g, h = (v + dt * (-drag / m - telemetry.G0 * math.sin(g)),
                       g + dt * (-telemetry.G0 * math.cos(g) / v),
                       h + dt * v * math.sin(g))
            t += dt
            if h < 1500 or v < 30:
                break


class Reports(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="caso_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_flight_cd_recovers_constant_cd(self):
        path = os.path.join(self.tmp, "vuelo.csv")
        synthetic_flight(path)
        s = {"file": path, "col_time": "0", "col_alt": "5", "col_speed": "11", "filter_col": "2",
             "filter_value": "TRACK", "mass_kg": "20", "aref_m2": "0.0179", "t_start": "", "smooth": ""}
        rows, note = telemetry.flight_cd(s, self.tmp)
        inner = sorted(r["cd_vuelo"] for r in rows[2:-2])
        self.assertAlmostEqual(inner[len(inner) // 2], 0.45, delta=0.45 * 0.02)
        self.assertIn("puntos", note)

    def test_flight_cd_from_accelerometer(self):
        path = os.path.join(self.tmp, "vuelo.csv")
        synthetic_flight(path)
        s = {"file": path, "col_time": "0", "col_alt": "5", "col_speed": "11", "col_accel": "7",
             "filter_col": "2", "filter_value": "TRACK", "mass_kg": "20", "aref_m2": "0.0179"}
        rows, note = telemetry.flight_cd(s, self.tmp)
        self.assertIn("acelerómetro", note)
        self.assertTrue(all(abs(r["cd_vuelo"] - 0.45) < 0.45 * 1e-3 for r in rows))

    def test_plan_paired_columns(self):
        caso = caso_de_prueba(self.tmp, """
            [caso]
            solver = openfoam
            id = M{mach}_{solver_tag}
            [sweep]
            mach,alt = 0.8,3000 1.8,1200
            solver_tag = a b
            """)
        rows = cfg.plan(caso)
        self.assertEqual([(r["mach"], r["alt"], r["solver_tag"]) for r in rows],
                         [("0.8", "3000", "a"), ("0.8", "3000", "b"), ("1.8", "1200", "a"), ("1.8", "1200", "b")])
        with self.assertRaises(ValueError):
            cfg.resolve(caso, dict(rows[0], alt="PENDIENTE"), machine(self.tmp))

    def test_missing_telemetry_config_is_a_note(self):
        rows, note = telemetry.flight_cd({"file": "x.csv"}, self.tmp)
        self.assertEqual(rows, [])
        self.assertIn("falta", note)

    def test_cost_curve_reference(self):
        caso = caso_de_prueba(self.tmp, """
            [caso]
            solver = nektar
            [report]
            kind = cost_curve
            group = mach
            family = solver stab
            reference = stab=dg
            """)
        rows = []
        for level in (0, 1, 2, 3):
            cd = 1.0 + 0.5 * (2.0 ** -level) ** 2
            rows.append({"id": "L%d" % level, "estado": "ok", "matriz.mach": "0.3", "matriz.stab": "dg",
                         "matriz.level": str(level), "matriz.P": "3", "corrida.solver": "nektar",
                         "resultado.cd": str(cd), "costo.horas_nucleo": str(4.0 ** level)})
        notes = []
        os.makedirs(caso.results_dir)
        report.cost_curve(caso, rows, caso.results_dir, notes)
        with open(os.path.join(caso.results_dir, "referencia.csv")) as f:
            ref = list(csv.DictReader(f))
        self.assertAlmostEqual(float(ref[0]["referencia"]), 1.0, 9)
        self.assertAlmostEqual(float(ref[0]["orden_observado"]), 2.0, 9)
        with open(os.path.join(caso.results_dir, "curva.csv")) as f:
            curve = list(csv.DictReader(f))
        self.assertAlmostEqual(float(curve[0]["error"]), 0.5, 9)

    def test_error_criteria(self):
        caso = caso_de_prueba(self.tmp, """
            [caso]
            solver = nektar
            [report]
            kind = errores
            case = caso
            columns = caso stab P h
            [criteria]
            V2 = linf <= 1e-12
            V3 = igual stab=dg tol=1e-10
            V4 = orden h variable=rhou margen=0.3
            """)

        def row(caso, stab, P, h, l2, linf):
            return {"id": "%s_%s_%s_%s" % (caso, stab, P, h), "estado": "ok", "matriz.caso": caso,
                    "matriz.stab": stab, "matriz.P": str(P), "matriz.h": str(h),
                    "errores.L2_rhou": str(l2), "errores.Linf_rhou": str(linf)}
        rows = [row("V2", "dg", 1, 1, 1e-14, 1e-13), row("V2", "svv", 1, 1, 1e-10, 1e-9),
                row("V3", "dg", 2, 1, 0.5, 0.5), row("V3", "svv", 2, 1, 0.5, 0.5),
                row("V4", "dg", 2, 0.2, 8e-3, 1), row("V4", "dg", 2, 0.1, 1e-3, 1),
                row("V4", "svv", 2, 0.2, 8e-3, 1)]
        os.makedirs(caso.results_dir)
        notes = []
        report.error_table(caso, rows, caso.results_dir, notes)
        with open(os.path.join(caso.results_dir, "errores.csv")) as f:
            got = {(r["caso"], r["stab"], r["h"]): r["cumple"] for r in csv.DictReader(f)}
        self.assertEqual(got[("V2", "dg", "1")], "si")
        self.assertEqual(got[("V2", "svv", "1")], "no")
        self.assertEqual(got[("V3", "dg", "1")], "")
        self.assertEqual(got[("V3", "svv", "1")], "si")
        self.assertEqual(got[("V4", "dg", "0.1")], "si")
        self.assertEqual(got[("V4", "svv", "0.2")], "pendiente")
        self.assertTrue(any("NO CUMPLE" in n for n in notes))

    def test_cd_mach_with_telemetry(self):
        path = os.path.join(self.tmp, "vuelo.csv")
        synthetic_flight(path)
        caso = caso_de_prueba(self.tmp, """
            [caso]
            solver = openfoam
            [report]
            kind = cd_mach
            [telemetry]
            file = %s
            filter_col = 2
            filter_value = TRACK
            col_time = 0
            col_alt = 5
            col_speed = 11
            mass_kg = 20
            aref_m2 = 0.0179
            """ % path)
        rows = [{"id": "m%s" % m, "estado": "ok", "matriz.mach": str(m), "corrida.solver": "openfoam",
                 "resultado.cd": "0.45"} for m in (0.2, 0.4, 0.8, 1.2, 1.8)]
        notes = []
        os.makedirs(caso.results_dir)
        report.cd_mach(caso, rows, caso.results_dir, notes)
        with open(os.path.join(caso.results_dir, "telemetria_cd.csv")) as f:
            lines = f.read().splitlines()
        self.assertGreater(len(lines), 10)
        self.assertFalse(any("menos de 2" in n for n in notes))


if __name__ == "__main__":
    unittest.main()
