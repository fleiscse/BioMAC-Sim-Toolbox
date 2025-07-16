% ======================================================================
%> @file @Gait2d_osim/Gait2d_osim.m
%> @brief Matlab class describing the Gait2d_osim model. This is based on
%> the Gait3d model and OpenSim's gait10dof18musc
%>
%> @author Ton, Eva, Anne, Marlies, Markus
%> @date May, 2022
% ======================================================================

%======================================================================
%> @brief The class describes the Gait2d_osim model
%>
%> @details
%> - The model is exactly OpenSim's gait10dof18musc model
%> - The code is copied from gait3d entirely, with one dimension removed
%======================================================================
classdef Gait2d_osim_tread < Gait2d_osim

    properties
        grf_part_left = zeros(1*(101-1),1);
        grf_part_right = zeros(1*(101-1),1);
        diff = zeros(1*(101-1),1);
        pd_part = zeros(1*(101-1),1);
    end


    properties (SetAccess = protected)
        %> Struct: Information (name, file, modified, sha256) on opensim model
       grf_delay = 4
       
       %for 1.2
       Kfx = 3.668972119324633e-06;%%%4.737902497799167e-06 %3.2389369371638603e-06%
       Kfy = -7.312253684547162e-07;%%%-9.399347670871730e-07 %-6.824653580500594e-07
       %Kgrf =   0.65008021190902970e-05 * 0.5;%3.224531316217335e-06%5.7275e-06%

       Kp =  0.588983902563203;%%%0.608387439669403%0.592267017990304%
   %3.804707484807214%2.2421%2.25%0.10131%4.25
       Kd = 0%;-0.007906228565342;

%%%%-0.009957138800960*0.7 %-0.006620853054847907
      
%   -0.014515311156096%-0.02523%-0.03164%-0.3034889%
       %Kpd =    0.121612336668859*0.5;%0.263247%0.2275%0.023701%0.2275


 %for 1.8
%         Kfy = -0.1722
%        Kgrf = 9.99e-06*0.6
% 
%        Kp = 3.1604
%        Kd = -0.0143
%        Kpd = 0.1619

       c = 0.01
    end
  
  

    methods

    

        %======================================================================
        %> @brief Function computing implicit differential equation for 3D musculoskeletal model
        %>
        %> @details
        %> This function calls the mex file of gait2d_osim.c:
        %> [f, dfdx, dfdxdot, dfdumus, dfdMextra]  = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);
        %>
        %> with the neural excitation umus and the extra torques Mextra.
        %>
        %> The dynamic residuals will be between fmin and fmax when inputs
        %> satisfy system dynamics : fmin <= f(x,dx/dt,umus,Mextra) <= fmax
        %>
        %> The last four outputs are optional and some computation time is saved if you do
        %> not request all of them.
        %>
        %> @param   obj         Gait2d_osim class object
        %> @param   x           Double array: State of the model (Gait2d_osim.nStates x 1)
        %> @param   xdot        Double array: State derivatives (Gait2d_osim.nStates x 1)
        %> @param   u           Double array: Controls of the model (Gait2d_osim.nControls x 1)
        %> @param   vBeltLeft   Double: Current speed of the left treadmill belt
        %> @param   vBeltRight  Double: Current speed of the right treadmill belt
        %>
        %> @retval  f       Double array: Dynamic residuals (Gait2d_osim.nConstraints x 1)
        %> @retval	dfdx	(optional) Double matrix: Transpose of Jacobian matrix df/dx 		(Gait2d_osim.nStates x Gait2d_osim.nConstraints)
        %> @retval	dfdxdot	(optional) Double matrix: Transpose of Jacobian matrix df/dxdot 	(Gait2d_osim.nStates x Gait2d_osim.nConstraints)
        %> @retval	dfdu	(optional) Double matrix: Transpose of Jacobian matrix df/du 		(Gait2d_osim.nControls x Gait2d_osim.nConstraints)
        %======================================================================
        function [f, dfdx, dfdxdot, dfdu] = getDynamics(obj,x,xdot,u, vBeltLeft, vBeltRight)
          %  stack = dbstack('-completenames');
           % if numel(stack) >= 2
            %    fprintf('Caller is %s line %d in file %s.\n', stack(2).name, stack(2).line, stack(2).file)
           % else
           %     fprintf('No caller.\n')
           % end
           
            % Get neural excitation
            idxu = obj.extractControl('u');
            umus = u(idxu);

            % Get extra moments (arm torques)
            Mextra = zeros(obj.nDofs, 1);
            idxTorque = obj.extractControl('torque');
            Mextra(obj.hidxTorqueDof) = obj.mExtraScaleFactor * u(idxTorque); % Scale them by obj.mExtraScaleFactor and assume that order is consistent.

            
            % Apply treadmill speed
            xdot(obj.idxcxCPleft) = xdot(obj.idxcxCPleft) + vBeltLeft; % add the current treadmill speed to the contaact point (instead of target speed)
            xdot(obj.idxcxCPright) = xdot(obj.idxcxCPright) + vBeltRight;

            % Get dynamics
            if nargout > 3

                [f, dfdx, dfdxdot, dfdumus, dfdMextra] = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);

                % Get dfdu from dfdumus and dfdMextra
                dfdu = zeros(obj.nControls, obj.nConstraints);
                dfdu(idxu, :) = dfdumus;
                dfdu(idxTorque, :) = dfdMextra(obj.hidxTorqueDof, :) * obj.mExtraScaleFactor;  % scaling has to be considered

            elseif nargout > 1
                [f, dfdx, dfdxdot] = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);
            else
                f = obj.hdlMEX('Dynamics',x,xdot,umus,Mextra);
            end

        end
        

    end

  

end



      