export output_pyna, generate_pyna_dirs
"""


generates a pyNA case in its preferred location. The case comprises a 
case directory (with subdirectories as shown in `generate_pyna_dirs`)
and the necessary data to run a noise assessment.

"""
function output_pyna(ac::TASOPT.aircraft=TASOPT.load_default_model())

    #check if ac has everything we need

    

    #generate the directory

    #save the data

    #hope for the best

end


"""


generates a directory structure at the specified pyna_path*"/cases" with 
for a specific case_name. considers an overwrite flag.

pyNA/cases
└── tasopt_default
    ├── aircraft
    ├── engine
    ├── output
    ├── shielding
    └── trajectory

"""
function generate_pyna_dirs(pyna_path, case_name="tasopt_default",
    overwrite = false)

    #check if base pyna path exists
    if !isdir(pyna_path)
        throw(ArgumentError("Base pyNA path could not be found: $pyna_path"))
    end

    # Create the full path to the case directory
    case_path = joinpath(pyna_path, "cases", case_name)

    # Check if the case directory already exists
    if isdir(case_path)
        if overwrite
            # If overwrite is true, delete the existing directory
            rm(case_path; force=true, recursive=true)
        else
            # If overwrite is false, throw an error
            throw(ArgumentError("Case directory $case_name already exists, and overwrite specified false."))
        end
    end

    # Create the case directory
    mkdir(case_path)

    # Create subdirectories within the case directory
    subdirectories = ["aircraft", "engine", "output", "shielding", "trajectory"]
    for subdir in subdirectories
        mkdir(joinpath(case_path, subdir))
    end
end
