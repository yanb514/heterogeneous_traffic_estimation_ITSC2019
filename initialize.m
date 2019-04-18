function U = initialize(model_param,test)
name = model_param.model_name;
switch name
    case 'True model'
        switch test
            case 1
                model_param.rm1 = 1.5; model_param.rm2 =1.0;
                model_param.vm1 = 1.8; model_param.vm2 = model_param.vm1;
            case 2
                model_param.rm1 = 1.5; model_param.rm2 = 1.0; % it was 1.8 and 1.0 satisfy eqn(7)
                model_param.vm1 = 1.8; model_param.vm2 = model_param.vm1;
        end
        
    case 'Creeping model'
        switch test
            case 1
                model_param.rm1 = 1.6; model_param.rm2 = 0.9;
                model_param.vm1 = 1.9; model_param.vm2 = model_param.vm1;
            case 2
                model_param.rm1 = 1.4; model_param.rm2 = 0.9; % it was 1.8 and 1.0 satisfy eqn(7)
                model_param.vm1 = 1.7; model_param.vm2 = model_param.vm1;
        end
        
    case 'LWR'
        switch test
            case 1
                model_param.rm1 = 1.8; model_param.rm2 = 1.0; %1.8 and 1.0 satisfy eqn(7)
                model_param.vm1 = 1.8; model_param.vm2 = model_param.vm1;
            case 2
                model_param.rm1 = 1.8; model_param.rm2 = 1.0; %1.8 and 1.0 satisfy eqn(7)
                model_param.vm1 = 1.8; model_param.vm2 = model_param.vm1;
        end
end
model_param.tfinal = 4;    % 150
model_param.len = 2;    % 50
model_param.x = linspace(0,model_param.len,40);
model_param.dx = model_param.x(2)-model_param.x(1);
model_param.dt = (model_param.x(2)-model_param.x(1))/2.16;
model_param.lambda = model_param.dt/model_param.dx;
model_param.t = 0:model_param.dt:model_param.tfinal;
model_param.M = length(model_param.t);
model_param.N = length(model_param.x);

%==========================================
switch test
    case 1   % overtaking
        thres1 = 2;
        thres2 = 0.6;
        thres1p5 = 0.9;
        thres3 = 1.4;
        freq = 0.07;
        switch name
            case 'True model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres1p5*model_param.N/model_param.len)) = .5;
                U(1,ceil(thres1p5*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:ceil(thres3*model_param.N/model_param.len)) = .7;
                U(2,ceil(thres3*model_param.N/model_param.len)+1:model_param.N) = eps;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.03+0.03)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.03+0.03)';
                model_param.d1r = 0+0*model_param.t;
                model_param.d2r = 0+0*model_param.t;
                
            case 'Creeping model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .6;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:ceil(thres3*model_param.N/model_param.len)) = .6;
                U(2,ceil(thres3*model_param.N/model_param.len)+1:model_param.N) = eps;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.09+0.09)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.09+0.09)';
                model_param.d1r = 0.05+0*model_param.t;
                model_param.d2r = 0.05+0*model_param.t;
                
            case 'LWR'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .5;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:ceil(thres3*model_param.N/model_param.len)) = .7;
                U(2,ceil(thres3*model_param.N/model_param.len)+1:model_param.N) = eps;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.00+0.00)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.00+0.00)';
                model_param.d1r = U(1,end)+0*model_param.t;
                model_param.d2r = U(2,end)+0*model_param.t;
        end
        
    case 2 % creeping
        thres1 = 2;
        thres2 = 0.7;
        freq = 0.07;
        switch name
            case 'True model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .5;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = .7;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.08+0.08)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.07+0.07)';
                model_param.d1r = 0.0+0*model_param.t;
                model_param.d2r = 1+0*model_param.t;
                
            case 'Creeping model'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .6;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = .6;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.03+0.03)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.03+0.03)';
                model_param.d1r = 0.1+0*model_param.t;
                model_param.d2r = 0.8+0*model_param.t;
                
            case 'LWR'
                U(1,1:thres1) = eps;
                U(1,thres1+1:ceil(thres2*model_param.N/model_param.len)) = .5;
                U(1,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = eps;
                U(2,1:ceil(thres2*model_param.N/model_param.len)) = eps;
                U(2,ceil(thres2*model_param.N/model_param.len)+1:model_param.N) = .7;
                
                model_param.d1l = (square(freq*[1:model_param.M])*0.00+0.00)';
                model_param.d2l = (square(freq*[1:model_param.M])*0.00+0.00)';
                model_param.d1r = 0+0*model_param.t;
                model_param.d2r = 0.999+0*model_param.t;
                
        end
end
end