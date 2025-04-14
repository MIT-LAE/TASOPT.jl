function compare_strings(user_input::String, reference::String)
    return lowercase(strip(user_input)) == lowercase(strip(reference))
end
