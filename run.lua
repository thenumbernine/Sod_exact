#!/usr/bin/env luajit
--[[
using:
1977 Sod "A Survey of Several Finite Difference Methods for Systems of Nonlinear Hyperbolic Conservation Laws"
http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/Chapter_6.pdf
https://en.wikipedia.org/wiki/Sod_shock_tube
--]]
require 'ext'
local matrix = require 'matrix'
local gnuplot = require 'gnuplot'

local n = 2000
local I = function(i) return i end
local C = function(c) return function() return c end end
local U = function(f) return range(n):map(f) end
local data = function(t) return table.map(t, U) end



--[[ idk why they use such eccentric values in the example
-- maybe they are closer to real-world values?
local min = -60
local max = 60
local t = 5000
local rhoL = 1e+5
local rhoR = 1.25e4
--]]
-- [[
local min = -5
local max = 5
local t = 1
local rhoL = 1
local rhoR = 1/8
--]]
local PL = 1
local PR = .1
local vL = 0
local vR = 0
local gamma = 5/3	--7/5
local muSq = (gamma - 1)/(gamma + 1)
local K = PL / rhoL^gamma

local i = matrix(range(n))
--local x = i:map(function(i) return (i-.5)/n * (max - min) + min end)
local x = (I - .5) / n * (max - min) + min

local CsL = math.sqrt(gamma * PL / rhoL)
local CsR = math.sqrt(gamma * PR / rhoR)

function solveP3()
	-- [[ Newton descent
	local symmath = require 'symmath'
	local P3 = symmath.var'P3'
	--[=[ Wikipedia 
	local v2 = (PL^((gamma-1)/2) - P3^((gamma-1)/2)) * symmath.sqrt( (1 - muSq*muSq) * PL^(1/gamma) / (muSq*muSq * rhoL) )
	local v4 = (P3 - PR) * symmath.sqrt( (1 - muSq) / (rhoR * (P3 + muSq * PR) ) )
	local f = v2 - v4
	--]=]
	-- [=[ http://www.itam.nsc.ru/flowlib/SRC/sod.f -- which is the same as Dullemon suggests in his text -- though sod.f uses a bisection intersection ...
	local f = -2*CsL*(1 - (P3/PL)^((-1 + gamma)/(2*gamma)))/(CsR*(-1 + gamma)) + (-1 + P3/PR)*((1 - muSq)/(gamma*(muSq + P3/PR)))^.5
	--]=]
	local df_dP3 = f:diff(P3)()	
	local f_func = f:compile{P3}
	local df_dP3_func = df_dP3:compile{P3}
	local P3 = .5 * (PL + PR)
	local epsilon = 1e-16	-- this is the limit for the sod.f before it oscillates
	while true do
		local dP3 = -f_func(P3) / df_dP3_func(P3)
		print(P3, dP3)
		if math.abs(dP3) <= epsilon then break end
		if not math.isfinite(dP3) then error('delta is not finite! '..tostring(dP3)) end
		P3 = P3 + dP3 
	end
	return P3
	--]]
end

local P3 = solveP3()
local P4 = P3

local rho3 = rhoL * (P3 / PL) ^ (1 / gamma)

-- Wikipedia:
--local v3 = vR + (P4 - PR) / math.sqrt( rhoR / 2 * ((gamma+1) * P3 + (gamma-1) * PR))
-- (although Wikipedia's root finding for P3 uses Dullemon's v3)
-- Dullemon:
--local v3 = (P4 - PR) * math.sqrt((1 - muSq) / (rhoR * P4 * muSq * PR))
-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
local v3 = vR + 2 * CsL / (gamma - 1) * (1 - (P3 / PL)^((gamma - 1)/(2*gamma)))
local v4 = v3

local rho4 = rhoR * (P4 + muSq * PR) / (PR + muSq * P4)

local vshock = v4 * rho4 / (rho4 - rhoR)
local vtail = CsL - v4 / (1 - muSq)

local v2 = function(x) return (1 - muSq) * (x/t + CsL) end
-- Dullemon:
--local rho2 = function(x) return (rhoL^gamma / (gamma * PL) * (v2(x) - x/t)^2)^(1/(gamma-1)) end
-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
local rho2 = function(x) return rhoL * (-muSq * (x / (CsL * t)) + (1 - muSq))^(2/(gamma-1)) end

-- Dullemon:
--local P2 = K * rho2^gamma
-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
local P2 = function(x) return PL * (-muSq * (x / (CsL * t)) + (1 - muSq)) ^ (2*gamma/(gamma-1)) end

local function makefunc(t)
	return table(t):map(function(ti,i)
		return type(ti) == 'function'and t[i]or C(t[i])
	end)
end

local rhos = makefunc{rhoL, rho2, rho3, rho4, rhoR}
local vs = makefunc{vL, v2, v3, v4, vR}
local Ps = makefunc{PL, P2, P3, P4, PR}

-- between regions 1 and 2
local s1 = -CsL	

-- between regions 2 and 3
-- book doesn't do this for me:
-- v2(s2) = v3(s2)
-- (1 - mu^2) (x/t + CsL) = v3
--  x/t + CsL = v3 / (1 - mu^2)
--  x/t = -CsL + v3 / (1 - mu^2)
--local s2 = s1 + v3 / (1 - mu^2)
-- http://www3.nd.edu/~gtryggva/CFD-Course/2013-Lecture-14.pdf says
-- says s2 = v4 - Cs3 ... so no (1 / mu^2) 
-- and it also doesn't give an equation for Cs3 ... so I'm assuming it's sqrt(gamma P3 / rho3)
--local Cs3 = math.sqrt(gamma * P3 / rho3)
--local s2 = -Cs3 + v4 
-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
local s2 = -vtail

local s3 = v3	-- between regions 3 and 4

-- between regions 4 and 5 ...
-- Wikipedia says this is Cs5
-- but is this necessarily the wave interface? 
-- it's certainly not in front of s3 ...
--local s4 = CsR
-- Dullemon 6.23
-- but this isn't just behind s3, it's behind everything other than s1
--local s4 = v4 * rho4 / (rho4 - rhoR)	
-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
local s4 = vshock

print('wavespeeds:',s1,s2,s3,s4)

local region = function(i)
	local x = x(i)
	local xi = x / t
	if xi < s1 then
		return 1
	elseif xi < s2 then
		return 2
	elseif xi < s3 then
		return 3
	elseif xi < s4 then
		return 4
	else
		return 5
	end
end

local rho = function(i) return rhos[region(i)](x(i)) end
local v = function(i) return vs[region(i)](x(i)) end
local P = function(i) return Ps[region(i)](x(i)) end
local eInt = P / (gamma - 1) 

local f = assert(io.open('sod_exact.txt', 'w'))
f:write'#x rho u P e\n'
local cols = {x, rho, v, P, eInt}
for i=1,n do
	for _,col in ipairs(cols) do
		f:write('\t', ('%.16f'):format(col(i)))
	end
	f:write'\n'
end
f:close()

local tmax = 5
local tres = 200
local ts = range(tres):map(function(i) return i/tres * tmax end)

-- [[ separate
for _,graph in ipairs{{rho=rho}, {v=v}, {P=P}} do
	local name,f = next(graph)
	gnuplot{
		output = name..'.png',
		style = 'data lines',
		data = data{x, f},
		{using='1:2', title=name},
	}
end
--]]
-- [[ together
gnuplot{
	output = 'results.png',
	style = 'data lines',
	data = data{x, rho, v, P, region},
	{using='1:2', title='rho'},
	{using='1:3', title='v'},
	{using='1:4', title='P'},
	{using='1:(($5-1)/4)', title='region'},
}
--]]
