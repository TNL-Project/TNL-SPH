def emptyTest( case_dir ):
    return 0, 0

configurations = []

# BIConsistent, MVT, Verlet, cs 25, 50, 100, DT filtering = 1, PST=Simple, p0 = 0
# BIConsistent, MVT, Verlet, cs 25, 50, 100, DT filtering = 5, PST=Simple, p0 = 0
# for n in [1, 5]:
#     for cs in [25, 50, 100]:
#         configurations.append({
#             "case-tag": f"spheric06-BIConsistent-MVT-Verlet-PST-DensityFilter-n{n}-cs{cs}",
#             "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#             "speed-of-sound": cs,
#             "density-filter": "ShepardFilter",
#             "filtering-steps-interval": n,
#             "evaluation-function": emptyTest
#         })

# # BIConsistent, MVT, Verlet, cs 25, 50, 100, DT filtering = 0, PST=Simple, TIC, p0 = 0
# for cs in [25, 50, 100]:
#     configurations.append({
#         "case-tag": f"spheric06-BIConsistent-MVT-Verlet-PST-TIC-cs{cs}",
#         "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#         "speed-of-sound": cs,
#         "grad-p" : "TIC",
#         "evaluation-function": emptyTest
#     })
# 
# # BIConsistent, MVT, Verlet, cs 25, 50, 100, DT filtering = 5, PST=Simple, TIC, p0 = 0
# for cs in [25, 50, 100]:
#     configurations.append({
#         "case-tag": f"spheric06-BIConsistent-MVT-Verlet-PST-TIC-DensityFilter-n5-cs{cs}",
#         "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#         "speed-of-sound": cs,
#         "density-filter": "ShepardFilter",
#         "filtering-steps-interval": 5,
#         "grad-p" : "TIC",
#         "evaluation-function": emptyTest
#     })

# # BIConsistent, MVT, symplectic, cs 25, 50, 100, 150, PST = None, p0 = 3
# for cs in [25, 50, 100, 150]:
#     configurations.append({
#         "case-tag": f"spheric06-BIConsistent-MVT-symplectic-p0_3-cs{cs}",
#         "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#         "time-integration": "SymplecticVerletScheme",
#         "speed-of-sound": cs,
#         "pst": "None",
#         "background-pressure": 3,
#         "evaluation-function": emptyTest
#     })

# # BIConsistent, MVT, Verlet, cs 25, 50, 100, 150, DT filtering = 1, PST= NONE, p0=3
# # BIConsistent, MVT, Verlet, cs 25, 50, 100, 150, DT filtering = 5, PST= NONE, p0=3
# for n in [1, 5]:
#     for cs in [25, 50, 100, 150]:
#         configurations.append({
#             "case-tag": f"spheric06-BIConsistent-MVT-symplectic-p0_3-DensityFilter-n{n}-cs{cs}",
#             "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#             "time-integration": "SymplecticVerletScheme",
#             "speed-of-sound": cs,
#             "density-filter": "ShepardFilter",
#             "filtering-steps-interval": n,
#             "pst": "None",
#             "background-pressure": 3,
#             "evaluation-function": emptyTest
#         })

# # BIConsistent, MGVT, Verlet, cs 25, 50, 100, DT,  PST=Simple, p0 = 0
# for cs in [25, 50, 100]:
#     configurations.append({
#         "case-tag": f"spheric06-BIConsistent-MGVT-Verlet-PST-cs{cs}",
#         "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#         "speed-of-sound": cs,
#         "viscous-term" : "PhysicalViscosity_MGVT",
#         "evaluation-function": emptyTest
#     })

# # BIConsistent, MGVT with MVT F eval, Verlet, cs 25, 50, 100, DT,  PST=Simple, p0 = 0
# for cs in [25, 50, 100]:
#     configurations.append({
#         "case-tag": f"spheric06-BIConsistent-MGVT-MVT_F_eval-Verlet-PST-cs{cs}",
#         "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#         "speed-of-sound": cs,
#         "viscous-term" : "PhysicalViscosity_MGVT",
#         "evaluation-function": emptyTest
#     })

# # BIConsistent, MVT and MGVT, Verlet, cs = 50, p0 = 0 (no features)
# for vt in ['MVT', 'MGVT']:
#     configurations.append({
#         "case-tag": f"spheric06-BIConsistent-{vt}-Verlet-cs50",
#         "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#         "speed-of-sound": 50,
#         "pst": "None",
#         "viscous-term" : f"PhysicalViscosity_{vt}",
#         "evaluation-function": emptyTest
#     })
# 
# # CONGERGENCE STUDY: BIConsistent, MVT, Verlet, cs 50 DT filtering = 5, PST=Simple, p0 = 0, dp 0.04 0.02 0.0125, 0.0075
# for dp in [0.04, 0.02, 0.0125, 0.0075]:
#     configurations.append({
#         "case-tag": f"spheric06-dp{dp}-BIConsistent-MVT-Verlet-PST-TIC-DensityFilter-n5-cs50",
#         "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
#         "dp": dp,
#         "speed-of-sound": 50,
#         "density-filter": "ShepardFilter",
#         "filtering-steps-interval": 5,
#         "evaluation-function": emptyTest
#     })

# CONVERGENCE STUDY: BIConsistent, MVT and MGVT, Verlet, cs = 50, p0 = 0 (no features)
for dp in [0.04, 0.02, 0.0125, 0.0075]:
    for vt in ['MVT', 'MGVT']:
        configurations.append({
            "case-tag": f"spheric06-dp{dp}-BIConsistent-{vt}-Verlet-cs50",
            "case": "WCSPH-BI/movingSquare2D_WCSPH-BI",
            "dp": dp,
            "speed-of-sound": 50,
            "pst": "None",
            "viscous-term" : f"PhysicalViscosity_{vt}",
            "evaluation-function": emptyTest
        })


from pprint import pprint
pprint(configurations)


