"""
     ___________________________________
    |                                   |
    |            fort2julia             |
    |              Authors:             |  
    |            Sicheng He             |
    |         Prakash Prashanth         |
    |            08/22/2022             |
    |___________________________________|
    Notes:
    The code is used to convert your .f77 file to .jl.
    It should give a good start point for further 
    editing. Further editing is required.

    Users need to do by themselves:
    1. Need to handle allocation, IO.
    2. "go to" stuff. 

    TODOs:
    1. The group operation is kind-of hacky. 
    May exist cleaner ways to handle that...
    2. rewrite "fix_patt_newline_math"
    3. 

"""


import re

# ==============================================
# Patterns
# ==============================================

# replace comment characters
patt_comment_c = re.compile(r"^c", flags=re.MULTILINE | re.IGNORECASE)
patt_comment_i = re.compile(r"!")

# Replace subroutines with functions
patt_sub = re.compile(r"subroutine")
# print(patt_sub.findall("      subroutine gasfun(igas,t, s,s_t, h,h_t, cp,r)"))

# Replace & line continuations with just a continuation
patt_cont = re.compile(r"(\s*)(&)(\s*)")

# Move any operator to previous line
patt_newline_math = re.compile(r"^(\s*)(&)(\s*)([(\-+*/)])")

# ** -> ^
patt_power = re.compile(r"\*\*")

# Trailing then
patt_then = re.compile(r"\s+then", flags=re.MULTILINE)

# endif enddo
patt_enddo = re.compile(r"enddo")
patt_endif = re.compile(r"endif")

# Replace .gt. .lt. etc
patt_gt = re.compile(r"\.gt\.")
patt_lt = re.compile(r"\.lt\.")
patt_eq = re.compile(r"\.eq\.")
patt_le = re.compile(r"\.le\.")
patt_ge = re.compile(r"\.ge\.")
patt_true = re.compile(r"\.true\.")
patt_false = re.compile(r"\.false\.")

# do loop
patt_do = re.compile(r"(do)(\s*)([0-9a-zA-Z_]*)(\s*)(=)(\s*)([0-9a-zA-Z_]*)(\s*)(,)(\s*)([0-9a-zA-Z_]*)")

# Call
# Assumption:
#           1. Function written in one line, e.g. call f1(x, y, \n & z), will not work
#           2. Max one function per line, e.g., call f1(x, z) call f2(x, y) \n, will not work.
#           3. No array showing up in the same line, e.g., call f1(x, z(j)), will not work.
patt_call = re.compile(r"(call)(\s+)([0-9a-zA-Z_]+)(\()([0-9a-zA-Z_,]+)(\))")

# Parenthesis to bracket for array
# Assumption:
#           1. Array written in one line, e.g. arr1(i, j, \n & k), will not work.
#           2. No function showing up in the same line, e.g., call f1(x, z(j)), will not work.
# Known issuse:
#           1. The allocation is messed up.
patt_arr_parenthesis = re.compile(r"([0-9a-zA-Z_]+)(\()([0-9a-zA-Z_\-+*/]+)(\))")

# Special fortran type: data
# 1: array handling
patt_data_bgn = re.compile(r"^(\s*)(data)(\s*)([0-9a-zA-Z_]+)(\s*)(/)")
patt_data_end = re.compile(r"(/)$")

# TODO []
# 2. (multiple) scalars handling
patt_data_scalar = re.compile(r"^(\s*)(data)(\s+)(([0-9a-zA-Z_])(\s*)(,)(\s*)*)([0-9a-zA-Z_]\s*)(\s*)")


# ==============================================
# Pattern fixes
# ==============================================


def fix_patt_newline_math(line):

    m = patt_newline_math.search(line)

    old_string = ""
    new_string = ""
    i = 0
    for substring in m.groups():
        old_string += substring
        if not i == 3:
            new_string += substring
        i += 1

    line = patt_newline_math.sub(new_string, line)

    return line


def fix_patt_cont(line):

    m = patt_cont.search(line)

    new_string = ""
    i = 0
    for substring in m.groups():
        if not i == 1:
            new_string += substring
        i += 1

    line = patt_cont.sub(new_string, line)

    return line


def fix_patt_comment_c(line):

    line = patt_comment_c.sub("#", line)

    return line


def fix_patt_comment_i(line):

    line = patt_comment_i.sub("#", line)

    return line


def fix_patt_sub(line):

    line = patt_sub.sub("function", line)

    return line


def fix_patt_power(line):

    line = patt_power.sub(r"^", line)

    return line


def fix_patt_then(line):

    line = patt_then.sub("", line)

    return line


def fix_patt_enddo(line):

    line = patt_enddo.sub("end", line)

    return line


def fix_patt_endif(line):

    line = patt_endif.sub("end", line)

    return line


def fix_patt_lt(line):

    line = patt_lt.sub("<", line)

    return line


def fix_patt_gt(line):

    line = patt_gt.sub(">", line)

    return line


def fix_patt_eq(line):

    line = patt_eq.sub("==", line)

    return line


def fix_patt_le(line):

    line = patt_le.sub("<=", line)

    return line


def fix_patt_ge(line):

    line = patt_ge.sub(">=", line)

    return line


def fix_patt_true(line):

    line = patt_true.sub("true", line)

    return line


def fix_patt_false(line):

    line = patt_false.sub("false", line)

    return line


def fix_patt_call(line):

    m = patt_call.search(line)

    new_string = ""
    for substring in m.groups():
        if not (substring == "call"):
            new_string += substring

    line = patt_call.sub(new_string, line)

    return line


def fix_patt_arr_parenthesis(line):

    if patt_call.search(line) == None:

        m_list = patt_arr_parenthesis.findall(line)

        for m in m_list:

            m_list = []
            for substring in m:
                m_list.append(substring)
            m_list[1] = "\("
            m_list[-1] = "\)"

            old_string = ""
            for substring in m_list:
                old_string += substring

            m_list[1] = "["
            m_list[-1] = "]"
            new_string = ""
            for substring in m_list:
                new_string += substring

            line = re.sub(old_string, new_string, line)

    return line


def fix_patt_do(line):

    m = patt_do.search(line)

    new_string = ""
    for substring in m.groups():
        if substring == "do":
            new_string += "for"
        elif substring == ",":
            new_string += ":"
        else:
            new_string += substring

    line = patt_do.sub(new_string, line)

    return line


def fix_data(lines):

    i_line = 0
    while 1:
        if i_line >= len(lines):
            break

        # Get the beginning and ending of the data
        # block.
        if patt_data_bgn.search(lines[i_line]) != None:

            i_line_data_bgn = i_line

            j_line = i_line + 1
            while 1:
                if j_line >= len(lines):
                    break

                if patt_data_end.search(lines[j_line]) != None:

                    i_line_data_end = j_line
                    break

                j_line += 1

            # Process it
            # Starting line
            m = patt_data_bgn.search(lines[i_line_data_bgn])

            new_string = ""
            for substring in m.groups():
                if substring == "data":
                    substring = ""
                elif substring == "/":
                    substring = "=["

                new_string += substring

            lines[i_line_data_bgn] = patt_data_bgn.sub(new_string, lines[i_line_data_bgn])

            # Ending line
            m = patt_data_end.search(lines[i_line_data_end])

            new_string = ""
            for substring in m.groups():
                if substring == "/":
                    substring = "]"

                new_string += substring

            lines[i_line_data_end] = patt_data_end.sub(new_string, lines[i_line_data_end])

        i_line = i_line + 1


# ==============================================
# Fix file and line
# ==============================================


def fix_file(input_file, output_file):
    try:
        f = open(input_file, "r")
    except IOError as err:
        print("{} could not be opened: {}".format(file, err))
        return 1

    # Read file to memory
    lines = f.readlines()

    # Features that detected in the line but affect previous fix here
    i_line = 0
    while 1:
        if i_line >= len(lines):
            break

        fix_previous_line(lines, i_line)
        i_line = i_line + 1

    # Fix multi-line data
    fix_data(lines)

    # Fix line-by-line
    i_line = 0
    while 1:
        if i_line >= len(lines):
            break

        lines[i_line] = fix_line(lines[i_line])
        i_line = i_line + 1

    f = open(output_file, "w")
    f.writelines(lines)
    f.close()


def fix_previous_line(lines, i_line):

    # Move any operator to previous line
    if patt_newline_math.search(lines[i_line]) != None:
        if patt_data_bgn.search(lines[i_line - 1]) == None:
            m = patt_newline_math.search(lines[i_line])
            operator_str = m.groups()[3]

            # Move the operator to the previous line
            lines[i_line - 1] = lines[i_line - 1][:-2] + operator_str + "\n"

            # clean this line
            lines[i_line] = fix_patt_newline_math(lines[i_line])


def fix_line(line):

    if patt_data_scalar.search(line) != None:

        print(line)

    # Fix the comment
    if patt_comment_c.search(line) != None:

        line = fix_patt_comment_c(line)

    if patt_comment_i.search(line) != None:

        line = fix_patt_comment_i(line)

    # Fix &
    if patt_cont.search(line) != None:

        line = fix_patt_cont(line)

    # Fix function
    if patt_sub.search(line) != None:

        line = fix_patt_sub(line)

    # Fix power
    if patt_power.search(line) != None:

        line = fix_patt_power(line)

    # Fix then
    if patt_then.search(line) != None:

        line = fix_patt_then(line)

    # Fix enddo endif
    if patt_enddo.search(line) != None:

        line = fix_patt_enddo(line)

    if patt_endif.search(line) != None:

        line = fix_patt_endif(line)

    # Fix logic
    if patt_lt.search(line) != None:

        line = fix_patt_lt(line)

    if patt_gt.search(line) != None:

        line = fix_patt_gt(line)

    if patt_eq.search(line) != None:

        line = fix_patt_eq(line)

    if patt_le.search(line) != None:

        line = fix_patt_le(line)

    if patt_ge.search(line) != None:

        line = fix_patt_ge(line)

    if patt_true.search(line) != None:

        line = fix_patt_true(line)

    if patt_false.search(line) != None:

        line = fix_patt_false(line)

    # Fix do
    if patt_do.search(line) != None:

        line = fix_patt_do(line)

    # Fix call
    if patt_call.search(line) != None:

        line = fix_patt_call(line)

    # Fix paranthesis
    if patt_arr_parenthesis.findall(line) != []:

        line = fix_patt_arr_parenthesis(line)

    return line


# fix_file("gasfun.f", "gasfun_fix.jl")
