function fractional_engine_weight!(ac)
    parg = ac.parg
    parg[igWeng] = parg[igWMTO] * parg[igfeng]
end