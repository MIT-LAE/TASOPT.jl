"""
    gaussn(nsiz, nn, z, r, nrhs)

Solves general nxn system in nn unknowns with arbitrary number (nrhs) of righthand sides. 
Assumes system is invertible. If it isn't, a divide by zero will result.
    Mark Drela  1984 

!!! details "🔃 Inputs and Outputs"
    **Input:**
    - `nsiz`: size?
    - `nn`: n * n number of unknowns in n x n system
    - `z`: coefficient matrix destroyed during solution process.
    - `r`: righthand side[s] replaced by the solution vector[s]
    - `nrhs`: Number of right hand sides

    **Output:**
    - `r`: solution of system

"""
function gaussn(nsiz, nn, z, r, nrhs)

  for np = 1:nn-1
    np1 = np + 1
    #
    #------ find max pivot index nx
    nx = np
    for n = np1:nn
      if ((abs(z[n, np]) - abs(z[nx, np])) > 0)
        nx = n
      end
    end
    #
    pivot = 1.0 / z[nx, np]
    #
    #------ switch pivots
    z[nx, np] = z[np, np]
    #
    #------ switch rows  normalize pivot row
    for l = np1:nn
      temp = z[nx, l] * pivot
      z[nx, l] = z[np, l]
      z[np, l] = temp
    end
    #
    for l = 1:nrhs
      temp = r[nx, l] * pivot
      r[nx, l] = r[np, l]
      r[np, l] = temp
    end
    #
    #------ forward eliminate everything
    for k = np1:nn
      ztmp = z[k, np]

      if !(ztmp == 0.0)
        for l = np1:nn
          z[k, l] = z[k, l] - ztmp * z[np, l]
        end
        for l = 1:nrhs
          r[k, l] = r[k, l] - ztmp * r[np, l]
        end
      end
    end

  end
  #
  #---- solve for last row
  for l = 1:nrhs
    r[nn, l] = r[nn, l] / z[nn, nn]
  end
  #
  #---- back substitute everything
  for np = nn-1:-1:1
    np1 = np + 1
    for l = 1:nrhs
      for k = np1:nn
        r[np, l] = r[np, l] - z[np, k] * r[k, l]
      end
    end
  end
  #
  return r
end # gauss