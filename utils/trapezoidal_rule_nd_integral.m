function out = trapezoidal_rule_nd_integral(x, mat, N)
%    N-Dimensional trapezoidal rule
%
%
% By: Mohammed S. Al-Rawi
% IEETA, University of Aveiro
% Portugal
% March, 8, 2013 
%
% Please send feedback if you have any issue with this function to
%  rawi707@yahoo.com
%
%
% This is a recursive function that computes an N-Dimensional integral 
% using the trapezoidal rule. 
% 
% ----- inputs:
% -mat is a N dimensional matrix that contains your data,
% -x is a cell array that contains the N vectors that
% correspond to the generation of mat
% - dim is the number of dimensions
%
%
% The best way to illustrate the use of this function is via an example,
% and here we go.
%
%
%
% 
% Example1: 3D data:
% x=0:.05:pi;
% y=0:.05:1;
% z=-1:.05:1;
% [xx,yy,zz] = ndgrid(x,y,z); 
% mat = yy.*sin(xx)+zz.*cos(xx);
% b{1}=x;b{2}=y; b{3}=z;
% out = trapezoidal_rule_nd_integral(b, mat, 3)
% 
% out =
% 
%     1.9987
%  Now, using Matlab's triple integration function
% F = @(x,y,z)y*sin(x)+z*cos(x);
% Q = triplequad(F,0,pi,0,1,-1,1)
% 
% Q =
% 
%     2.0000
%
%
% Example 2: 4D data
%
%
% x1=0:.01:1;
% x2=0:.05:1;
% x3=-1:.07:1;
% x4 = 0:.1:pi;
% [xx1,xx2,xx3, xx4] = ndgrid(x1,x2,x3,x4); 
% mat = xx1.*sin(xx2)+xx3.*cos(xx4);
% b{1}=x1;b{2}=x2; b{3}=x3;b{4}=x4;
% out = trapezoidal_rule_nd_integral(b, mat, 4)
% out = trapezoidal_rule_nd_integral(b, mat, 4)
%
% out =
%
%    1.3946
%
% Oh, am too lazy to check the answer, would you do that for me please,
% anyone out there?
%
% From https://de.mathworks.com/matlabcentral/fileexchange/40690-n-dimensional-trapezoidal-integral
% Comment by Eike Petersen, Aug 2021: I did check this somewhat thoroughly, and in all cases I looked at
% the results agreed with those of integral / integral2 / integral3. Of course, the accuracy will be reduced.
%
mat = trapz(x{N}, mat, N);
if N==1
    out=mat;
    return;
end
out = trapezoidal_rule_nd_integral(x, mat, N-1);
% N.B.  I am not sure how the error of this integration would grow with
% higher dimensions, please check any numerical analysis textbook, or any
% website that discusses this issue. The accuracy of trapz is less than
% other methods even in 1D. Overall, you need to be careful of the output
% you will get with this function, although it is very handy. 
