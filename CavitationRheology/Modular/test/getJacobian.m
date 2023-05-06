function J = getJacobian(name)
    if name == "DiffEqn1Compt"
        J = @JacobianDiffEqn1Compt;
    elseif name == "DiffConv1Compt"
        J = @JacobianDiffConv1Compt;
    elseif name == "Vdp"
        J = @JacobianVdp;
    end
end

function J = JacobianVdp(params, y, x, helpers)
    mu = params('mu');
    J = [0, 1; -2*mu*(y(1)*y(2))-1, mu*(1-y(1)^2)];
    %J = [0, -2*mu*(y(1)*y(2))-mu; 1, mu*(1-y(1)^2)];
end

function J = JacobianDiffEqn1Compt(params, y, x, helpers)
    D = params("D");
    numX = params("numX");
    dx = params("dx");
    %numX = (xSpan(2)-xSpan(1)) / dx;
    J = zeros(numX, numX);
%     for i = 1:numX
%         if i == 1
%             J(i,i) = -2*D / (dx^2);
%             J(i,i+1) = 2*D / (dx^2);
%         elseif (i > 1 && i < numX)
%             J(i,i-1) = D / (dx^2);
%             J(i,i) = -2*D / (dx^2);
%             J(i,i+1) = D / (dx^2);
%         elseif i == numX
%             J(i,i-1) = 2*D / (dx^2);
%             J(i,i) = -2*D / (dx^2);
%         end
%     end
    for i = 1:numX
        if i == 1
            J(i,i) = -2*D / (dx^2);
            J(i,i+1) = 2*D / (dx^2);

            J(i+1,i) = 2*D / (dx^2);
        elseif (i > 1 && i < numX)
            J(i,i-1) = D / (dx^2);
            J(i,i) = -2*D / (dx^2);
            J(i,i+1) = D / (dx^2);
            
            J(i-1,i) = D / (dx^2);
            J(i+1,i) = D / (dx^2);
        elseif i == numX
            J(i,i-1) = 2*D / (dx^2);
            J(i,i) = -2*D / (dx^2);

            J(i-1,i) = 2*D / (dx^2);
        end
    end

end

function J = JacobianDiffConv1Compt(params, y, x, helpers)
    %rSpan = params.rSpan;
    %dr = params.dr;
    D = params("D");
    r = helpers{1};
    v = helpers{2};
    dvdr = helpers{3};
    numR = (length(r));
    J = zeros(numR, numR);
%     for i = 1:numR
%         if i == 1
%             dr = r(2)-r(1);
%             J(i,i) = -6*D / (dr^2);
%             J(i,i+1) = 6*D / (dr^2);
%             
%             J(i+1,i) = 6*D / (dr^2);
%         elseif (i > 1 && i < numR)
%             dr = r(i+1)-r(i);
%             J(i,i-1) = (-2*D / (r(i) * dr^2)) + D/(dr^2) + v(i)/(dr);
%             J(i,i) = (-2*D / (r(i) * dr^2)) - 2*D / (dr^2) - 2*v(i)/(r(i)) - dvdr(i) - v(i)/dr;
%             J(i,i+1) = D / (dr^2);
% 
%             J(i-1,i) = (-2*D / (r(i) * dr^2)) + D/(dr^2) + v(i)/(dr);
%             J(i+1,i) = D / (dr^2);
%         elseif i == numR
%             dr = r(i)-r(i-1);
%             J(i,i-1) = 6*D / (dr^2);
%             J(i,i) = -6*D / (dr^2);
% 
%             J(i-1,i) = 6*D / (dr^2);
%         end
%     end
    for i = 1:numR
        if i == 1
            dr = r(2)-r(1);
            J(i,i) = -6*D / (dr^2);
            J(i,i+1) = 6*D / (dr^2);
        elseif (i > 1 && i < numR)
            dr = r(i+1)-r(i);
            J(i,i-1) = (-2*D / (r(i) * dr^2)) + D/(dr^2) + v(i)/(dr);
            J(i,i) = (-2*D / (r(i) * dr^2)) - 2*D / (dr^2) - 2*v(i)/(r(i)) - dvdr(i) - v(i)/dr;
            J(i,i+1) = D / (dr^2);
        elseif i == numR
            dr = r(i)-r(i-1);
            J(i,i-1) = 6*D / (dr^2);
            J(i,i) = -6*D / (dr^2);
        end
    end
end