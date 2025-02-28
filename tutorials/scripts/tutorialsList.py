# tutorials list for python validation script

import json

tutorials = [
    {
        "solver": "anisoImpesFoam",
        "cases": [
            { "case": "injectionAniso_case1", "type": "sample",
              "result": "postProcessing/sampleDict/2500/acrossFlow_Sb.csv" },
            { "case": "injectionAniso_case2", "type": "sample",
              "result": "postProcessing/sampleDict/2500/acrossFlow_Sb.csv" }
        ]
    },
    {
        "solver": "impesFoam",
        "cases": [
            { "case": "Buckley-Leverett/BrooksAndCorey", "type": "sample",
              "result": "postProcessing/sampleDict/20000/acrossFlow_Sb.csv" },
            { "case": "Buckley-Leverett/VanGenuchten", "type": "sample",
              "result": "postProcessing/sampleDict/30000/acrossFlow_Sb.csv" },
            { "case": "capillarityValidation/BrooksAndCorey", "type": "sample",
              "result": "postProcessing/sampleDict/1000/acrossFlow_Sb.csv" },
            { "case": "capillarityValidation/VanGenuchten", "type": "sample",
              "result": "postProcessing/sampleDict/1000/acrossFlow_Sb.csv" },
            { "case": "injectionExtraction/injection", "type": "sample",
              "result": "postProcessing/sampleDict/25000/acrossFlow_Sb.csv" },
            { "case": "injectionExtraction/extraction", "type": "sample", 
              "result": "postProcessing/sampleDict/25000/acrossFlow_Sb.csv" }
        ]
    },
    {
        "solver": "darcyFoam",
        "cases": [
            { "case": "SPE10",  "type": None, "result": None }
        ]
    },
    {
        "solver": "groundwater2DFoam",
        "cases": [
            { "case": "steadyFlow", "type": "steady", "label": "potential",
              "result": "potential_probes.dat"},
            { "case": "transientFlow", "type": "probe", "label": "potential",
              "result": "postProcessing/probes/0/potential", "col": 5}
        ]
    },
    {
        "solver": "groundwaterFoam",
        "cases": [
            { "case": "1Dinfiltration", "type": "sample",
              "result": "postProcessing/sampleDict/50000/acrossFlow_h.csv"},
            { "case": "1Dinfiltration_dualPorosity", "type": "sample",
              "result": "postProcessing/sampleDict/3456/acrossFlow_theta_thetaMatrix.csv"},
            { "case": "steadyFlow", "type": "sample",
              "result": "postProcessing/sampleDict/159/acrossFlow_theta.csv"},
            { "case": "transientFlow", "type": "transient", "label": "outflow",
              "result": "waterMassBalance.csv", "col": 2}
        ]
    },
    {
        "solver": "groundwaterTransport2DFoam",
        "cases": [
            { "case": "coupled", "type": "probe", "label": "concentration",
              "result": "postProcessing/probes/0/C", "col": 5}
        ]
    },
    {
        "solver": "groundwaterTransportFoam",
        "cases": [
            { "case": "1Dinfiltration_dualPorosity", "col": 5, "type": "sample",
              "result": "postProcessing/sampleDict/3456/acrossFlow_C_CMatrix_theta_thetaMatrix.csv"},
            { "case": "coupled", "type": "transient", "label": "Coutflow",
              "result": "CmassBalance.csv", "col": 3}
        ]
    },
    {
        "solver": "porousScalarTransport2DFoam",
        "cases": [
            { "case": "transport", "type": "probe", "label": "C_outflow",
             "result": "postProcessing/probes/0/C", "col": 5 },
        ]
    },
    {
        "solver": "porousScalarTransportFoam",
        "cases": [
            { "case": "1DeventFlux_Euler", "type": "transient", "label": "Coutflow",
              "result": "CmassBalance.csv", "col": 3},
            { "case": "1DeventFlux_backward", "type": "transient", "label": "Coutflow",
              "result": "CmassBalance.csv", "col": 3},
            { "case": "1DeventFlux_CrankNicolson", "type": "transient", "label": "Coutflow",
              "result": "CmassBalance.csv", "col": 3},
            { "case": "1DeventFlux_multispecies", "type": "transient", "label": "Coutflow",
              "result": "C1massBalance.csv", "col": 3},
            { "case": "transport", "type": "transient", "label": "Coutflow",
              "result": "CmassBalance.csv", "col": 3}
        ]
    }
]

# Convertit la structure en une chaîne JSON formatée
json_output = json.dumps(tutorials, indent=4)

