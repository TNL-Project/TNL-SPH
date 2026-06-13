CONFIGURATIONS = {

    # Test simple refined block
    "dummy-center": {
        "description": "Centered refinement in unit square",
        "dp": 0.002,
        "h_coef": 1.75,
        "box_x": 1.0,
        "box_y": 1.0,
        "fine_factor": 0.5,
        "fine_x_min": 0.25,
        "fine_x_max": 0.75,
        "fine_y_min": 0.25,
        "fine_y_max": 0.75,
    },

    # testcase - dambreak
    "dambreak2D-corner": {
        "description": "Dam break with refinement in bottom-right corner",
        "dp": 0.002,
        "h_coef": 1.75,
        "box_x": 1.61,
        "box_y": 0.6,
        "fine_factor": 0.5,
        "fine_x_min": 1.2,
        "fine_x_max": None,
        "fine_y_min": None,
        "fine_y_max": 0.3,
    },

    # Test corners
    # (note: the subdomain limits are randomly selected)
    "dummy-pp-corner": {
        "description": "Centered refinement in unit square",
        "dp": 0.02,
        "h_coef": 1.75,
        "box_x": 1.0,
        "box_y": 1.0,
        "fine_factor": 0.5,
        "fine_x_min": 0.63,
        "fine_x_max": None,
        "fine_y_min": 0.85,
        "fine_y_max": None,
    },

    "dummy-pp-corner-h2": {
        "description": "Centered refinement in unit square",
        "dp": 0.02,
        "h_coef": 2,
        "box_x": 1.0,
        "box_y": 1.0,
        "fine_factor": 0.5,
        "fine_x_min": 0.65,
        "fine_x_max": None,
        "fine_y_min": 0.85,
        "fine_y_max": None,
    },

    "dummy-pm-corner": {
        "description": "Centered refinement in unit square",
        "dp": 0.02,
        "h_coef": 1.75,
        "box_x": 1.0,
        "box_y": 1.0,
        "fine_factor": 0.5,
        "fine_x_min": 0.75,
        "fine_x_max": None,
        "fine_y_min": None,
        "fine_y_max": 0.75,
    },

    "dummy-mm-corner": {
        "description": "Centered refinement in unit square",
        "dp": 0.02,
        "h_coef": 1.75,
        "box_x": 1.0,
        "box_y": 1.0,
        "fine_factor": 0.5,
        "fine_x_min": None,
        "fine_x_max": 0.34,
        "fine_y_min": None,
        "fine_y_max": 0.27,
    },

    "dummy-mp-corner": {
        "description": "Centered refinement in unit square",
        "dp": 0.02,
        "h_coef": 1.75,
        "box_x": 1.0,
        "box_y": 1.0,
        "fine_factor": 0.5,
        "fine_x_min": None,
        "fine_x_max": 0.25,
        "fine_y_min": 0.75,
        "fine_y_max": None,
    },

    "dummy-mp-corner-h2": {
        "description": "Centered refinement in unit square",
        "dp": 0.02,
        "h_coef": 2,
        "box_x": 1.0,
        "box_y": 1.0,
        "fine_factor": 0.5,
        "fine_x_min": None,
        "fine_x_max": 0.25,
        "fine_y_min": 0.75,
        "fine_y_max": None,
    },

}
