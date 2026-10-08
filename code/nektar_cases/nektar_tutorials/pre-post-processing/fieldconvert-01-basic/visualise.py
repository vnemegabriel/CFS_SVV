import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
import matplotlib

from NekPy.LibUtilities import SessionReader, ShapeType
from NekPy.SpatialDomains import MeshGraph
from NekPy.MultiRegions import ExpList

def triangulate_2d(session_files, fieldfile, varName, subdivs=6):
    session = SessionReader.CreateInstance(['tmp'] + session_files)
    graph = MeshGraph.Read(session)
    graph.SetExpansionInfosToEvenlySpacedPoints(subdivs)
    exp = ExpList(session, graph)
    exp.LoadField(fieldfile, varName)

    x, y = exp.GetCoords()
    tris = []

    for i in range(exp.GetExpSize()):
        elmt = exp.GetExp(i)
        nq0 = elmt.GetBasis(0).GetNumPoints()
        nq1 = elmt.GetBasis(1).GetNumPoints()
        offset = exp.GetPhys_Offset(i)

        for j in range(nq0-1):
            for k in range(nq1-1):
                tris.append([offset + nq0 * k + j,
                             offset + nq0 * k + j+1,
                             offset + nq0 * (k+1) + j])
                tris.append([offset + nq0 * k + j+1,
                             offset + nq0 * (k+1) + j+1,
                             offset + nq0 * (k+1) + j])

    return matplotlib.tri.Triangulation(x, y, tris), exp

def plot_mesh(exp, linewidth=1, linestyle='solid', color='k'):
    x, y = exp.GetCoords()
    lines = []
    for i in range(exp.GetExpSize()): 
        elmt = exp.GetExp(i)
        nq0 = elmt.GetBasis(0).GetNumPoints()
        nq1 = elmt.GetBasis(1).GetNumPoints()
        offset = exp.GetPhys_Offset(i)

        for j in range(nq0-1):
            lines.append([
                (x[offset + j], y[offset + j]), (x[offset + j + 1], y[offset + j + 1])
            ])
        for j in range(nq1-1):
            lines.append([
                (x[offset + j*nq0], y[offset + j*nq0]), (x[offset + (j+1)*nq0], y[offset + (j+1)*nq0])
            ])
            tmp = offset + nq0-1
            lines.append([
                (x[tmp + j*nq0], y[tmp + j*nq0]), (x[tmp + (j+1)*nq0], y[tmp + (j+1)*nq0])
            ])

        if elmt.DetShapeType() == ShapeType.Triangle:
            continue

        for j in range(nq0-1):
            tmp = offset + nq0 * (nq1-1)
            lines.append([
                (x[tmp + j], y[tmp + j]), (x[tmp + j + 1], y[tmp + j + 1])
            ])

    segments = LineCollection(lines, colors=color, linestyles=linestyle, linewidths=linewidth)
    ax = plt.gca()
    ax.add_collection(segments)

#tri, exp = triangulate_2d(['adr-geometry.xml', 'adr-conditions.xml'], 'adr-geometry.fld', 'u')
#plt.tricontourf(tri, exp.GetPhys())
#plot_mesh(exp)
#plt.show()
