#!/usr/bin/lua
-- ----------------------------------------------------------------- --
--      This Lua5 script is Copyright (c) 2010, Peter J Billam       --
--                        www.pjb.com.au                             --
--                                                                   --
--   This script is free software; you can redistribute it and/or    --
--          modify it under the same terms as Lua5 itself.           --
-- ----------------------------------------------------------------- --
local Version = '1.0  for Lua5'
local VersionDate  = '31jul2010';
local Synopsis = [[
program_name [options] [filenames]
]]
local iarg=1; while arg[iarg] ~= nil do
	if string.sub(arg[iarg],1,1) ~= "-" then break end
	local first_letter = string.sub(arg[iarg],2,2)
	if first_letter == 'v' then
		local n = string.gsub(arg[0],"^.*/","",1)
		print(n.." version "..Version.."  "..VersionDate)
		os.exit(0)
	elseif first_letter == 'c' then
		whatever()
	else
		local n = string.gsub(arg[0],"^.*/","",1)
		print(n.." version "..Version.."  "..VersionDate.."\n\n"..Synopsis)
		os.exit(0)
	end
	iarg = iarg+1
end

local RK = require 'RungeKutta'
-- use Test::Simple tests => 6;

local func_evals = 0
local i_test     = 0
local n_failed   = 0
local n_passed   = 0
function dydt(t, y)
	local dydt_f = {}
	dydt_f[1] = y[2]
	dydt_f[2] = 0.0 - y[1]
	func_evals = func_evals + 1
	return dydt_f
end

twopi = 2.0 * 3.141592653589
algorithms = {RK.rk2, RK.rk4, RK.rk4_classical, RK.rk4_ralston, epsilon, errors}
passmark0 = { 0.2,   .0004,  .0015,               .0015, .0001, .0003 }
passmark1 = { 0.04,  .00004, .0006,               .0006, .00001, .00001 }
for i=1,4 do
	algorithm = algorithms[i]
	i_test = i_test + 1
	n = 16
	dt= twopi / n
	y = {0,1}; t=0

	for j=1, n do t, y = algorithm(y, dydt, t, dt) end
	local err0 = math.abs(y[1]);
	local err1 = math.abs(y[2]-1.0)
	if err0 < passmark0[i] and err1 < passmark1[i] then print('passed '..i)
	else print('failed '..i_test)
	end
end

algorithm = RK.rk4_auto
local t_midpoint; local y_midpoint
for k, mode in ipairs({'epsilon','errors'}) do
	local skip_to_next = false
	i_test = i_test + 1
	y = {0,1}; t = 0; dt = 0.1
	local i = 0
	local epsilon
	if mode == 'epsilon' then epsilon = .0001
	else errors = {.01, .0001}; epsilon = errors;
	end
	func_evals = 0
	while t+dt < twopi do
		i = i + 1
		t, dt, y = RK.rk4_auto(y, dydt, t, dt, epsilon )
-- print('t = '..t..'  dt = '..dt..'  y = '..y[1]..' '..y[2])
		t_midpoint, y_midpoint = RK.rk4_auto_midpoint()
		if (func_evals > 500) then
			print('rk4_auto('..mode..') failed, '..func_evals..' func evals')
			skip_to_next = true ; break
		end
	end
	if not skip_to_next then
		i = i + 1; dt = twopi-t;
		t, y = RK.rk4(y, dydt, t, dt)
-- print('y = '..y[1]..' '..y[2])
		local err0 = math.abs(y[1]);
		local err1 = math.abs(y[2]-1.0)
		if err0 < passmark0[i_test] and err1 < passmark1[i_test] then
			print('passed rk4_auto('..mode..')')
		else
			print('failed rk4_auto('..mode..')')
		end
	end
end

--[[
__END__

=pod

=head1 NAME

test.lua - Lua script to test RungeKutta.lua

=head1 SYNOPSIS

 perl test.pl

=head1 DESCRIPTION

This script tests Math::RungeKutta.pm

=head1 AUTHOR

Peter J Billam http://www.pjb.com.au/comp/contact.html

=head1 SEE ALSO

Math::RungeKutta.pm , http://www.pjb.com.au/ , perl(1).

=cut

]]
