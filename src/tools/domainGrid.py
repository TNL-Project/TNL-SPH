from contextlib import redirect_stdout

def write_domain_grid(setup: dict, filename: str) -> None:
    search_radius = setup["search_radius"]
    grid_size_x = round(setup["domain_size_x"] / search_radius)
    grid_size_y = round(setup["domain_size_y"] / search_radius)
    grid_size_z = round(setup.get("domain_size_z", 0) / search_radius)

    grid_origin_x = setup["domain_origin_x"]
    grid_origin_y = setup["domain_origin_y"]
    grid_origin_z = setup.get("domain_origin_z", 0)

    n_cells = grid_size_x * grid_size_y * max(grid_size_z, 1)

    with open(filename, "w") as f:
        with redirect_stdout(f):
            print("# vtk DataFile Version 3.0")
            print("vtk output")
            print("ASCII")
            print("DATASET STRUCTURED_POINTS")
            print("DIMENSIONS ", grid_size_x + 1, " ", grid_size_y + 1, " ", grid_size_z + 1)
            print("ASPECT_RATIO ", search_radius, " ", search_radius, " ", search_radius)
            print("ORIGIN ", grid_origin_x, " ", grid_origin_y, " ", grid_origin_z)
            print("CELL_DATA ", n_cells)
            print("SCALARS GridSector int 1 ")
            print("LOOKUP_TABLE default")
            for _ in range(n_cells):
                print(0)
