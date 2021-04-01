import pandas as pd


data = []

v0 = 5.98
v0sink = 11.96

for ivol_sink in [-0.5, -0.25, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5]:
    for ivol_source in [-0.5, -0.25, -0.2, -0.1, -0.05, -0.025, 0, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5]:
        for ipower in [50, 60, 70, 80, 90, 100]:

            hp.sink_in.vdot = v0 * (ivol_sink + 1)
            hp.sink_out.vdot = hp.sink_in.vdot
            hp.source_in.vdot = v0source * (ivol_source + 1)
            hp.source_out.vdot = hp.source_in.vdot

            data.append([ivol_sink, ivol_source, ipower, v0 * (ivol_sink + 1), v0source * (ivol_source + 1) ])

nosol = pd.DataFrame(data, columns=["sink Vdot", "Source vdot", "Ipower"])

with pd.ExcelWriter("runtest.xlsx") as writer:
    nosol.to_excel(writer, sheet_name="NoSolution")
