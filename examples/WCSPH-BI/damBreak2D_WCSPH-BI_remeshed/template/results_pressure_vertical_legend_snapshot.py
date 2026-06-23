# state file generated using paraview version 6.0.1
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.Set(
    ViewSize=[1536, 795],
    InteractionMode='2D',
    AxesGrid='Grid Axes 3D Actor',
    OrientationAxesVisibility=0,
    CenterOfRotation=[0.8050000071525574, 0.30000001192092896, 0.0],
    CameraPosition=[0.8824475714355292, 0.29781872637491635, 3.3192450147143995],
    CameraFocalPoint=[0.8824475714355292, 0.29781872637491635, 0.0],
    CameraFocalDisk=1.0,
    CameraParallelScale=0.5867658118773873,
    UseColorPaletteForBackground=0,
    Background=[1.0, 1.0, 1.0],
    OSPRayMaterialLibrary=materialLibrary1,
)

# init the 'Grid Axes 3D Actor' selected for 'AxesGrid'
renderView1.AxesGrid.Set(
    Visibility=1,
    XTitle='',
    YTitle='',
    ZTitle='',
    FacesToRender=7,
    GridColor=[0.7803921699523926, 0.7803921699523926, 0.7803921699523926],
    ShowGrid=1,
    AxesToLabel=7,
    XLabelFontFamily='Times',
    XLabelFontSize=36,
    YLabelFontFamily='Times',
    YLabelFontSize=36,
    ZLabelFontFamily='Times',
    ZLabelFontSize=36,
    XAxisNotation='Fixed',
    XAxisPrecision=1,
    XAxisUseCustomLabels=1,
    XAxisLabels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
    YAxisNotation='Fixed',
    YAxisPrecision=1,
    YAxisUseCustomLabels=1,
    YAxisLabels=[0.0, 0.2, 0.4, 0.6],
)

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1536, 795)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
fluid_0000000_particlesvtk = LegacyVTKReader(registrationName='fluid_0.000000_particles.vtk', FileNames=['/home/tomas/Work/sph-projects/fresh/tnl-sph/examples/WCSPH-BI/damBreak2D_WCSPH-BI/results/fluid_0.000000_particles.vtk'])

# create a new 'Legacy VTK Reader'
boundary_0000000_particlesvtk = LegacyVTKReader(registrationName='boundary_0.000000_particles.vtk', FileNames=['/home/tomas/Work/sph-projects/fresh/tnl-sph/examples/WCSPH-BI/damBreak2D_WCSPH-BI/results/boundary_0.000000_particles.vtk'])

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from fluid_0000000_particlesvtk
fluid_0000000_particlesvtkDisplay = Show(fluid_0000000_particlesvtk, renderView1, 'GeometryRepresentation')

# get 2D transfer function for 'Pressure'
pressureTF2D = GetTransferFunction2D('Pressure')
pressureTF2D.Set(
    ScalarRangeInitialized=1,
    Range=[0.0, 3000.0, 0.0, 1.0],
)
fluid_0000000_particlesvtkDisplay.RescaleTransferFunctionToDataRange(True, True)

# get color transfer function/color map for 'Pressure'
pressureLUT = GetColorTransferFunction('Pressure')
pressureLUT.Set(
    TransferFunction2D=pressureTF2D,
    RGBPoints=[
        # scalar, red, green, blue
        0.0, 0.0, 0.0, 0.5625,
        333.3330000000002, 0.0, 0.0, 1.0,
        1095.2385, 0.0, 1.0, 1.0,
        1476.1905, 0.5, 1.0, 0.5,
        1857.1425, 1.0, 1.0, 0.0,
        2619.0480000000002, 1.0, 0.0, 0.0,
        3000.0, 0.5, 0.0, 0.0,
    ],
    ColorSpace='RGB',
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
fluid_0000000_particlesvtkDisplay.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Pressure'],
    LookupTable=pressureLUT,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
fluid_0000000_particlesvtkDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
fluid_0000000_particlesvtkDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
fluid_0000000_particlesvtkDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show data from boundary_0000000_particlesvtk
boundary_0000000_particlesvtkDisplay = Show(boundary_0000000_particlesvtk, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
boundary_0000000_particlesvtkDisplay.Set(
    Representation='Outline',
    ColorArrayName=['POINTS', ''],
    DiffuseColor=[0.0, 0.0, 0.0],
    LineWidth=2.0,
    DataAxesGrid='Grid Axes Representation',
    PolarAxes='Polar Axes Representation',
)

# init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
boundary_0000000_particlesvtkDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
boundary_0000000_particlesvtkDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
boundary_0000000_particlesvtkDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for pressureLUT in view renderView1
pressureLUTColorBar = GetScalarBar(pressureLUT, renderView1)
pressureLUTColorBar.Set(
    WindowLocation='Any Location',
    Position=[0.857121772670025, 0.24523073998384093],
    Title='Pressure',
    ComponentTitle='',
    TitleFontFamily='Times',
    TitleFontSize=20,
    LabelFontFamily='Times',
    LabelFontSize=20,
    ScalarBarThickness=20,
    ScalarBarLength=0.5124944948750856,
    DrawScalarBarOutline=1,
    ScalarBarOutlineColor=[0.0, 0.0, 0.0],
    ScalarBarOutlineThickness=2,
    AddRangeLabels=0,
)

# set color bar visibility
pressureLUTColorBar.Visibility = 1

# show color legend
fluid_0000000_particlesvtkDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Pressure'
pressurePWF = GetOpacityTransferFunction('Pressure')
pressurePWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 3000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=renderView1,
    Cues=timeAnimationCue1,
    AnimationTime=0.0,
)

# initialize the animation scene

# ----------------------------------------------------------------
# restore active source
SetActiveSource(None)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/nightly/python/
##--------------------------------------------
