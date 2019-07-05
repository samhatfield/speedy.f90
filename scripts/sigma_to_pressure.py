from argparse import ArgumentParser
from iris import load_cube, load
from iris.exceptions import ConstraintMismatchError
from iris import Constraint, save
from iris.coords import DimCoord
from cf_units import Unit
import numpy as np
from os.path import splitext

# Pressure levels to interpolate to
p_levels = np.array([25.0, 100.0, 200.0, 300.0, 500.0, 700.0, 850.0, 950.0], dtype=np.float32)


def interpolate_σ_to_p(p_levels, σ_pressure_levels, input):
    output = input.copy()
    for lat in range(input.shape[1]):
        for lon in range(input.shape[2]):
            output[:,lat,lon] \
                = np.interp(p_levels, σ_pressure_levels[:,lat,lon], input[:,lat,lon])

    return output


# Parse command line arguments
parser = ArgumentParser(description="Converts the given file from sigma level to pressure level")
parser.add_argument("filename", type=str, help="File name of input")
args = parser.parse_args()

# Load surface pressure cube
try:
    pₛ = load_cube(args.filename, "surface_air_pressure")
except ConstraintMismatchError:
    print(f"Couldn't find surface pressure cube in file {args.filename}")
    print("Is it named correctly? (surface_air_pressure)")
    raise SystemExit

# Load other cubes
cubes = load(args.filename, Constraint(cube_func=lambda n: n.long_name != "surface_air_pressure"))
print("Processing the following cubes...")
print(cubes)

# Get sigma coordinates
σ_levels = cubes[0].coord("atmosphere_sigma_coordinate").points

print(f"Processing {args.filename} from sigma levels " +
      f"\n{', '.join(str(σ) for σ in σ_levels)}" +
      f"\nto pressure levels\n{', '.join(str(p) + ' hPa' for p in p_levels)}")

# Construct pressure level Iris DimCoord
pressure_level_coord = DimCoord(p_levels, units=Unit("hPa"), long_name="air_pressure",
                                var_name="lev")

# Copy cubes and replace σ level DimCoord with pressure level DimCoord
for cube in cubes:
    cube.remove_coord("atmosphere_sigma_coordinate")
    cube.add_dim_coord(pressure_level_coord, 1)

for t, slice in enumerate(pₛ.slices_over("time")):
    print(f"Converting time slice {t+1} of {pₛ.shape[0]}")

    # Convert sigma levels to equivalent pressure levels (making sure to convert to hPa)
    σ_pressure_levels = slice.data/100.0 * σ_levels[:,None,None]

    for cube in cubes:
        cube.data[t,...] = interpolate_σ_to_p(p_levels, σ_pressure_levels, cube.data[t,...])

# Save cubes to file
new_filename = f"{splitext(args.filename)[0]}_p.nc"
print(f"Saving to {new_filename}")
save(cubes, new_filename)
