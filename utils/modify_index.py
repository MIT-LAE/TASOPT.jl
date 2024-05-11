import re
import click

@click.command()
@click.option('--input_file', '-i', required=True, 
               help="Path to index file that you want to de-fortranify")
@click.option('--variable_names', '-vars', required=True, type=str,
              help="Enter a list of index name to be removed separated by commas.")
@click.option('--output_file', '-o', default=None, 
              help="Give output file name")
def remove_variable_and_renumber_indices(input_file, output_file,
                                          variable_names):
    print("Input file:",input_file)
    if output_file == None:
        output_file = input_file + "_modified"
    print("Output file:",output_file)
    # print(variable_names)
    variable_names = variable_names.split(',')
    print("Indices to remove:",variable_names,"\n")

    # Read the contents of the file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    for variable_name in variable_names:
        var_start = variable_name[0:2] #Store this so can only act on vars of this category

        # Find the index of the variable to be removed
        removed_index = None
        removed_line_no = None
        n_lines = len(lines)
        for i, line in zip(range(0,n_lines), lines):
            match = re.search(rf'\b{variable_name}\b\s*=\s*(\d+)\s*', line)
            if match:
                # print(i,line)
                removed_index = int(match.group(1))
                removed_line_no = i
                print(f"Removed Ind {variable_name} [{removed_index}] from '{var_start}' category")
                lines[i] = "### Good Riddance ###\n"
                break

        # Renumber the indices for the remaining variables in that category
        for i, line in zip(range(removed_line_no,n_lines), lines[removed_line_no:-1]):
            # Only match variables of the same category
            match = re.search(rf'\b{var_start}\w+\s*=\s*(\d+)\s*', line)
            if match:
                index = int(match.group(1))
                if index > removed_index:
                    # print("Index",index)
                    new_index = index - 1
                    lines[i] = re.sub(rf'\b{index}\b', str(new_index), line)

    # Write the modified list to a new text file
    print("\nWriting to file...")
    with open( output_file, 'w') as file:
        file.write(''.join(lines))
    print("Done!")

if __name__ == '__main__':
    remove_variable_and_renumber_indices()