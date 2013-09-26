-- refinement parameters
min_level = 2
max_level = 4
output_prefix = "solution"

-- connectivity info
connectivity = "brick"
brick =
{
  nx = 5,
  ny = 5,
  periodic_x = 1,
  periodic_y = 1,
}

-- set up the domain
Lx = 25
Ly = 25
function connectivity_vertices(x, y, z)
  xout = Lx*x
  yout = Ly*y
  zout = 0
  return xout,yout,zout
end

function refinement_function(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  if level < min_level then
    return 1
  -- elseif level >= max_level or level > treeid/4 then
  --   return 0
  -- else
  --   return 1
  else
    return 0
  end
end

function element_order(
  x0,y0,z0,x1,y1,z1,
  x2,y2,z2,x3,y3,z3,
  level, treeid)

  -- N = treeid%3+1
  N = 3

  return N
end

-- initial values
rho = 1
lam = 1
mu  = 1
S11 = 0
S22 = 0
S33 = 0
S12 = 0
S13 = 0
S23 = 0

-- plane wave parameters
A_p = 1
p_p = {1,0,0}
d_p = p_p
c_p = math.sqrt((lam+2*mu)/rho);
k_p = 2*math.pi/(5*Lx);

A_s = 1
p_s = {1,0,0}
d_s = {0,1,1}
c_s = math.sqrt(mu/rho);
k_s = 2*math.pi/(5*Ly);

function v(x1,x2,x3,t,i)
  return -A_p*d_p[i]*c_p*k_p*math.cos(k_p*(p_p[1]*x1+p_p[2]*x2+p_p[2]*x3-c_p*t))
         -A_s*d_s[i]*c_s*k_s*math.cos(k_s*(p_s[1]*x1+p_s[2]*x2+p_s[2]*x3-c_s*t))
end

function v1(x,y,z,t)
  return v(x,y,z,t,1)
end
function v2(x,y,z,t)
  return v(x,y,z,t,2)
end
function v3(x,y,z,t)
  return v(x,y,z,t,3)
end

function S(x1,x2,x3,t,i,j)
  S_p = mu*(d_p[i]*p_p[j]+d_p[j]*p_p[i])
  S_s = mu*(d_s[i]*p_s[j]+d_s[j]*p_s[i])
  if i == j then
    S_p = S_p + lam*(d_p[1]*p_p[1]+d_p[2]*p_p[2]+d_p[3]*p_p[3])
    S_s = S_s + lam*(d_s[1]*p_s[1]+d_s[2]*p_s[2]+d_s[3]*p_s[3])
  end
  return  A_p*k_p*S_p*math.cos(k_p*(p_p[1]*x1+p_p[2]*x2+p_p[2]*x3-c_p*t))
         +A_s*k_s*S_s*math.cos(k_s*(p_s[1]*x1+p_s[2]*x2+p_s[2]*x3-c_s*t))
end

function S11(x,y,z,t)
  return S(x,y,z,t,1,1)
end
function S22(x,y,z,t)
  return S(x,y,z,t,2,2)
end
function S33(x,y,z,t)
  return S(x,y,z,t,3,3)
end
function S12(x,y,z,t)
  return S(x,y,z,t,1,2)
end
function S13(x,y,z,t)
  return S(x,y,z,t,1,3)
end
function S23(x,y,z,t)
  return S(x,y,z,t,2,3)
end



-- time stepper to use
lsrk_method  = "KC54"

function time_step_parameters(dt)
  tend  = Lx*5/c_s
  tout  = tend/100
  tdisp = tend/100

  dt      = 0.5*dt
  noutput = math.ceil(tout / dt)
  dt      = tout / noutput
  ndisp   = tdisp / dt
  nstesp  = tend / dt
  return dt,nstesp, ndisp, noutput
end
