
function compare(a, d)

      eps0 = 5.0e-9
      eps1 = 5.0e-5

      ss = "**"
      if (a == d || (abs(a) < eps0 && abs(d) < eps0))

            ss = "  "

      elseif (abs(a) > abs(d))

            if (abs(1.0 - d / a) < eps1)
                  ss = "  "
            end

      elseif (abs(d) > abs(a))
            if (abs(1.0 - a / d) < eps1)
                  ss = "  "
            end
      end

      return ss
end

