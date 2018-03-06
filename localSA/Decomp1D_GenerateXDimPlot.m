function [xdim_plot] = Decomp1D_GenerateXDimPlot( xdim_plot, n, dx )
%   Generate the x-dimension matrix for plotting results
%
%    Calling Sequence:
%    [xdim_plot] = Decomp1D_GenerateXDimPlot( xdim_plot, dx )
%
%    Input:   xdim_plot - initialized xdim_plot matrix from main code
%             dx        - matrix of individual element lengths
%
%    Output:  xdim_plot - x-dimension matrix for plotting results
%

for i=1:n
    if(i==1)
        xdim_plot(i) = 0;
    else   
        xdim_plot(i) = xdim_plot(i-1)+dx(i-1);
    end
end
    