using Revise
using TASOPT

if !@isdefined ac
    
    ac = load_default_model()
    size_aircraft!(ac)
end

pyna_path = "C:\\Users\\jnsgn\\workspace\\pyNA\\argonaut\\pyNA"
output_pyna(pyna_path, ac, overwrite = true)
