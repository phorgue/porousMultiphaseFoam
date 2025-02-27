# tutorials list for python validation script

import json

tutorials = {
    "tutorials": [
        {
            "solver": "anisoImpesFoam",
            "cases": [
                { "case": "injectionAniso_case1" },
                { "case": "injectionAniso_case2" }
            ],
            "dataTypes": [
                { "dataType": ".csv" },
                { "dataType": ".csv" },
            ]
        },
        {
            "solver": "impesFoam",
            "cases": [
                { "case": "Buckley-Leverett/BrooksAndCorey" },
                { "case": "Buckley-Leverett/VanGenuchten" },
                { "case": "capillarityValidation/BrooksAndCorey" },
                { "case": "capillarityValidation/VanGenuchten" },
                { "case": "injectionExtraction/injection" },
                { "case": "injectionExtraction/extraction" }
            ],
            "dataTypes": [
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" },
            ]
        },
        {
            "solver": "darcyFoam",
            "cases": [
                { "case": "SPE10" }
            ],
            "dataTypes": [
                { "dataType": ".csv" }
            ]
        },
        {
            "solver": "groundwater2DFoam",
            "cases": [
                { "case": "steadyFlow" },
                { "case": "transientFlow" }
            ],
            "dataTypes": [
                { "dataType": "potential" },
                { "dataType": "potential" }
            ]
        },
        {
            "solver": "groundwaterFoam",
            "cases": [
                { "case": "1Dinfiltration" },
                { "case": "1Dinfiltration_dualPorosity" },
                { "case": "steadyFlow" },
                { "case": "transientFlow" }
            ],
            "dataTypes": [
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" }
            ]
        },
        {
            "solver": "groundwaterTransport2DFoam",
            "cases": [
                { "case": "coupled" }
            ],
            "dataTypes": [
                { "dataType": "C" },
                { "dataType": "potential" }
            ]
        },
        {
            "solver": "groundwaterTransportFoam",
            "cases": [
                { "case": "1Dinfiltration_dualPorosity" },
                { "case": "coupled" }
            ],
            "dataTypes": [
                { "dataType": ".csv" },
                { "dataType": ".csv" }
            ]
        },
        {
            "solver": "porousScalarTransport2DFoam",
            "cases": [
                { "case": "transport" }
            ],
            "dataTypes": [
                { "dataType": "C" },
            ]
        },
        {
            "solver": "porousScalarTransportFoam",
            "cases": [
                { "case": "1DeventFlux_Euler" },
                { "case": "1DeventFlux_backward" },
                { "case": "1DeventFlux_CrankNicolson" },
                { "case": "1DeventFlux_multispecies" },
                { "case": "transport" }
            ],
            "dataTypes": [
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" },
                { "dataType": ".csv" }
            ]
        }
    ]
}

# Convertit la structure en une chaîne JSON formatée
json_output = json.dumps(tutorials, indent=4)

