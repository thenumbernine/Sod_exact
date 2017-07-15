#!/usr/bin/env luajit
-- from http://www.mpia.de/homes/dullemon/lectures/fluiddynamics08/chap_6_numhyd_riemann_1.pdf 
require 'ext'
local matrix = require 'matrix'
local gnuplot = require 'gnuplot'

local n = 2000

local min = -60
local max = 60
local x = matrix{n}:lambda(function(i)
	return (i-.5)/n * (max - min) + min
end)

--[[ idk why they use such eccentric values in the example
-- maybe they are closer to real-world values?
local t = 5000
local rhoL = 1e+5
local rhoR = 1.25e4
--]]
-- [[
local t = 1
local rhoL = 1
local rhoR = 1/8
--]]
local PL = 1
local PR = .1
local vL = 0
local vR = 0
local gamma = 7/5
local m = (gamma - 1)/(gamma + 1)
local K = PL / rhoL^gamma

local P4 = 0 --numerically computed
local P3 = P4

local CsL = math.sqrt(gamma * PL / rhoL)
local v2 = function(x) return (1 - m^2) * (x/t + CsL) end
local v3 = (P4 - PR) * math.sqrt( (1-m^2) / (rhoR * (P4 + m^2 * PR)) )
local v4 = v3

local rho4 = rhoR * (P4 * m^2 * PR) / (PR + m^2 * P4)
local rho3 = rhoL * (P3 / PL) ^ (1 / gamma)
local rho2 = function(x) return (rhoL^gamma / (gamma * PL) * (v2 - x/t)^2)^(1/(gamma-1)) end
local P2 = K * rho2^gamma

local rhos = {rhoL, rho2, rho3, rho4, rhoR}
local vs = {vL, v2, v3, v4, vR}
local Ps = {PL, P2, P3, P4, PR}

local s1 = -CsL * t
-- book doesn't do this for me:
-- v2(s2) = v3(s2)
-- (1 - m^2) (x/t + CsL) = v3
--  x/t + CsL = v3 / (1 - m^2)
--  x = -t*CsL + v3 / (1 - m^2)
--  x = s1 + v3 / (1 - m^2)
local s2 = s1 + v3 / (1 - m^2)
local s3 = v3
local s4 = v4

local i = x:map(function(x)
	if x < s1 then
		return 1
	elseif x < s2 then
		return 2
	elseif x < s3 then
		return 3
	elseif x < s4 then
		return 4
	else
		return 5
	end
end)

local rho = i:map(function(i) return rhos[i] end)
local v = i:map(function(i) return vs[i] end)
local P = i:map(function(i) return Ps[i] end)

--[[ separate
for _,graph in ipairs{{rho=rho}, {v=v}, {P=P}} do
	local name,data = next(graph)
	gnuplot{
		output = name..'.png',
		style = 'data lines',
		data = {x, data},
		{using='1:2', title=name},
	}
end
--]]
-- [[ together
gnuplot{
	output = 'results.png',
	style = 'data lines',
	data = {x, rho, v, P},
	{using='1:2', title='rho'},
	{using='1:3', title='v'},
	{using='1:4', title='P'},
}
--]]
