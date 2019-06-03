classdef ModeSolver
    properties
        lam_, r_, n1_, n2_
        lam, r, k0, n1, n2
        ndim
        indeps_, indeps_vals, indeps_keys
        indeps
        vvmax, HE11flag, n1lam, n2lam
        neffmlist_, neff
        modes, fields
        res
    end
    methods
        function obj = ModeSolver(lam, r, n1, n2, varargin)
            p = inputParser;
            addParameter(p, 'indeps', nan);
            addParameter(p, 'vvmax', 1);
            addParameter(p, 'HE11', true);
            addParameter(p, 'neffm', nan);
            addParameter(p, 'n1lam', nan);
            addParameter(p, 'n2lam', nan);
            addParameter(p, 'Dw', false);
            addParameter(p, 'neffmin', -12);
            parse(p, varargin{:});
            
            mkdir('tempFunctions')
            obj.indeps_ = p.Results.indeps;
            
            if ~iscell(obj.indeps_)
                obj.ndim = 0;
                [obj.indeps_keys, obj.indeps_vals] = deal(nan, nan);
            else
                obj.ndim = 1;
                [obj.indeps_keys, obj.indeps_vals] = deal(obj.indeps_{:});
            end
            
            obj.lam_ = lam;
            obj.r_ = r;
            obj.n1_ = n1;
            obj.n2_ = n2;
            
            obj.HE11flag = p.Results.HE11;
            obj.vvmax = p.Results.vvmax;
            obj.n1lam = p.Results.n1lam;
            obj.n2lam = p.Results.n2lam;
            if obj.vvmax > 1
                obj.HE11flag = false;
            end
            
            temp = linspace(0.1, 1, 101);
            neffm = p.Results.neffm;
            if isnan(neffm)
                obj.neffmlist_ = [10.^linspace(p.Results.neffmin, -1, 1001), temp(2:end)];
            else
                obj.neffmlist_ = neffm;
            end
            temp = obj.basicParaCal(obj.indeps_vals, 'neff', obj.neffmlist_);
            [temp2, obj.indeps, neffmlist] = deal(temp{:});
            [obj.lam, obj.k0, obj.r, obj.n1, obj.n2] = deal(temp2{:});
            obj.neff = obj.n2 + (obj.n1-obj.n2).*neffmlist;

            HEEHstrs = {'HE', 'EH'};
            obj.modes = struct(...
                'HEEH', struct(...
                    'name', @(v, i) [HEEHstrs{mod(i-1, 2)+1}, num2str(v), num2str(ceil(i/2))],...
                    'mode', @obj.eigenEq_HE_EH_vm,...
                    'dev', @obj.dev_eigenEq_HE_EH_vm...
                ),...
                'TE', struct(...
                    'name', @(i) ['TE0', num2str(i)],...
                    'mode', @obj.eigenEq_TE0m,...
                    'dev', @obj.dev_eigenEq_TE0m...
                ),...
                'TM', struct(...
                    'name', @(i) ['TM0', num2str(i)],...
                    'mode', @obj.eigenEq_TM0m,...
                    'dev', @obj.dev_eigenEq_TM0m...
                )...
            );
            obj.fields = obj.fieldDefine();

            obj.res = obj.dispersiveCurveCalc({obj.k0, obj.r, obj.n1, obj.n2});
            obj.res = obj.fieldCalc(obj.res);
            addpath(genpath(pwd));
            if p.Results.Dw
                obj = obj.DwCalc();
            end 
        end
        
        function delete(~)
            rmdir('tempFunctions', 's')
        end

        function res = dispersiveCurveCalc(obj, paras, varargin)
            p = inputParser;
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;
            
            if obj.HE11flag
                resMode = containers.Map({'HEEH'}, {{obj.calcFiberMode(paras, obj.neff, obj.modes.HEEH.mode, 'v', 1)}});
            else
                resMode = containers.Map({'HEEH', 'TE', 'TM'},...
                    {cellfun(@(vv) obj.calcFiberMode(paras, obj.neff, obj.modes.HEEH.mode, 'v', vv), num2cell(1:obj.vvmax), 'UniformOutput', false),...
                    obj.calcFiberMode(paras, obj.neff, obj.modes.TE.mode),...
                    obj.calcFiberMode(paras, obj.neff, obj.modes.TM.mode)},...
                    'UniformValues',false...
                );
            end

            res = struct();
            for key0 = keys(resMode)
                key = key0{1};
                vals = resMode(key);
                if strcmp(key, 'HEEH')
                    for v = 1:length(vals)
                        val = vals{v};
                        inters = obj.dispersiveInter(val, key, 'v', v, 'delam', delam);
                        for i = 1:length(inters)
                            res.(obj.modes.(key).name(v, i)) = struct('neff', inters{i});
                        end
                    end
                else
                    inters = obj.dispersiveInter(vals, key, 'delam', delam);
                    for i = 1:length(inters)
                        res.(obj.modes.(key).name(i)) = struct('neff', inters{i});
                    end
                end
            end
        end

        function paras = basicParaCal(obj, indeps, varargin)
            p = inputParser;
            addParameter(p, 'neff', 1);
            addParameter(p, 'usage', 'default');
            addParameter(p, 'output', 'paras');
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;
            neff = p.Results.neff;
            
            if strcmp(p.Results.usage, 'default')
                if obj.ndim == 1
                    [indeps, neff] = ndgrid(indeps, p.Results.neff);
                end
            end
            temp = {};
            paras = {obj.lam_, obj.r_, obj.n1_, obj.n2_};
            for i = 1:length(paras)
                para = paras{i};
                if isa(para, 'function_handle')
                    if iscell(indeps)
                        indeps = indeps{1};
                    end
                    temp{i} = para(indeps);
                else
                    temp{i} = para;
                end
            end
            [lam, r, n1, n2] = deal(temp{:}); %#ok<*PROPLC>
            lam = lam + delam;
            if ~isnumeric(obj.n1lam)
                n1 = obj.n1lam(lam);
            end
            if ~isnumeric(obj.n2lam)
                n2 = obj.n2lam(lam);
            end
            k0 = 2*pi./lam;
            if strcmp(p.Results.output, 'simple')
                paras = {lam, k0, r, n1, n2};
            else
                paras = {{lam, k0, r, n1, n2}, indeps, neff};
            end
        end

        function val = dispersiveInter(obj, val, modeName, varargin)
            p = inputParser;
            addParameter(p, 'v', nan);
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;
            
            temp = {@obj.dispersiveInter0D, @obj.dispersiveInter1D};
            getIsos = temp{obj.ndim+1};
            val = getIsos(val);

            for i=1:length(val)
                val{i} = sortrows(val{i});
            
                if obj.ndim > 0
                    vall = val{i};
                    [~, ia, ~] = unique(vall(:, 1:end-1), 'rows', 'stable');
                    val{i} = vall(ia, :);
                end

                [jj, ~] = size(val{i});
                for j=1:jj
                    if obj.ndim == 0
                        val{i}(j, end) = obj.findRoot(nan, val{i}(j, end), modeName, 'usage', 'findRoot', 'v', p.Results.v, 'delam', delam);
                    else
                        val{i}(j, end) = obj.findRoot(val{i}(j, 1:end-1), val{i}(j, end), modeName, 'usage', 'findRoot', 'v', p.Results.v, 'delam', delam);
                    end
                end

                vall = val{i};
                inan = ~isnan(vall(:, end));
                val{i} = vall(logical(inan'), :);
            end
            if obj.ndim > 0
                val = val(cellfun(@(x) length(x)>4, val));
            end

            for i=1:length(val)
                if obj.ndim == 0
                    val{i} = @() val{i};
                elseif obj.ndim == 1
                    vall = val{i}';
                    val{i} = @(x) interp1(vall(1, :), vall(2, :), x, 'spline', nan);
                end
            end
 
            if obj.ndim == 0
                [~, ind] = sort(cellfun(@(x) x(), val), 'descend');
            else
                [~, ind] = sort(cellfun(@(x) sum(x(obj.indeps_vals), 'omitnan'), val), 'descend');
            end
            val = val(ind);
            if obj.HE11flag
                val = val(1);
            end
        end

        function res = findRoot(obj, indeps, neff, modeName, varargin)
            p = inputParser;
            addParameter(p, 'usage', 'default');
            addParameter(p, 'delam', 0);
            addParameter(p, 'v', nan);
            parse(p, varargin{:});
            delam = p.Results.delam;
            
            mode = obj.modes.(modeName).mode;
            dev = obj.modes.(modeName).dev;
            temp = obj.basicParaCal(indeps, 'neff', neff, 'usage', p.Results.usage, 'output', 'simple', 'delam', delam);
            [~, k0, r, n1, n2] = deal(temp{:});
            fvalue = @(x) obj.calcFiberMode({k0, r, n1, n2}, x, mode, 'v', p.Results.v);
            prime = @(x) obj.calcDevFiberMode({k0, r, n1, n2}, x, dev, 'v', p.Results.v);
            func = @(x) deal(fvalue(x), prime(x));
            options = optimoptions('fsolve','Display','off','SpecifyObjectiveGradient',true, 'StepTolerance', 1e-12, 'FunctionTolerance', 1e-12, 'OptimalityTolerance', 1e-200);
            res = fsolve(func, neff, options);
            if abs(res-neff) > 0.2*(neff-1)
                res = nan;
            end
        end

        function res = dispersiveInter0D(obj, val)
            res = {};
            ii = 1;
            for i = 2:length(val)
                if sign(val(i-1))*sign(val(i))<1
                    res{ii} = obj.neff(i);
                    ii = ii + 1;
                end
            end
        end

        function csvs = dispersiveInter1D(obj, val)
            cs = contourc(obj.indeps_vals, obj.neff(1, :), val', [0 0]);
            csvs = {};
            ii = 1;
            i = 1;
            ns = size(cs);
            while i < ns(2)
                len = cs(2, i);
                i = i + 1;
                csvs{ii} = cs(:, i:i+len-1)';
                ii = ii + 1;
                i = i + len;
            end
        end

        function res = eigenEq_HE_EH_vm(~, ~, U, W, n1, n2, ~, v)
            besseljd = @(v, z)  1./2.*(besselj(v-1, z) - besselj(v+1, z));
            besselkd = @(v, z) -1./2.*(besselk(v-1, z) + besselk(v+1, z));
            res = (besseljd(v, U)./(U.*besselj(v, U))+besselkd(v, W)./(W.*besselk(v, W))) .* (besseljd(v, U)./(U.*besselj(v, U))+n2.^2./n1.^2.*besselkd(v, W)./(W.*besselk(v, W))) - v.^2 .* (1./U.^2+1./W.^2) .* (1./U.^2+n2.^2./n1.^2.*1./W.^2);
        end

        function res = dev_eigenEq_HE_EH_vm(~, ~, U, W, n1, n2, k0, neff, r, v)
            res = (k0.*neff.*r.*(-16.*n1.^2.*v.^2.*(U.^2 + W.^2).*besselj(v,U).^3.*besselk(v,W) - 16.*v.^2.*(n2.^2.*U.^2 + n1.^2.*W.^2).*besselj(v,U).^3.*besselk(v,W) + n1.^2.*U.^2.*W.*(2.*besselj(v,U).*(besselj(-1 + v,U) - besselj(1 + v,U)) + U.*(besselj(-1 + v,U) - besselj(1 + v,U)).^2 + U.*besselj(v,U).*(-besselj(-2 + v,U) + 2.*besselj(v,U) - besselj(2 + v,U))).*                 (W.*(besselj(-1 + v,U) - besselj(1 + v,U)).*besselk(v,W) - U.*besselj(v,U).*(besselk(-1 + v,W) + besselk(1 + v,W))) +                 U.^2.*W.*(2.*besselj(v,U).*(besselj(-1 + v,U) - besselj(1 + v,U)) + U.*(besselj(-1 + v,U) - besselj(1 + v,U)).^2 +                 U.*besselj(v,U).*(-besselj(-2 + v,U) + 2.*besselj(v,U) - besselj(2 + v,U))).*                 (n1.^2.*W.*(besselj(-1 + v,U) - besselj(1 + v,U)).*besselk(v,W) - n2.^2.*U.*besselj(v,U).*(besselk(-1 + v,W) + besselk(1 + v,W)))))./                 (8..*n1.^2.*sqrt(n1.^2 - neff.^2).*U.^5.*W.^2.*besselj(v,U).^3.*besselk(v,W)) +                 (k0.*neff.*r.*(16.*n2.^2.*v.^2.*(U.^2 + W.^2).*besselj(v,U).*besselk(v,W).^3 + 16.*v.^2.*(n2.^2.*U.^2 + n1.^2.*W.^2).*besselj(v,U).*besselk(v,W).^3 +                 n2.^2.*U.*W.^2.*(W.*(besselj(-1 + v,U) - besselj(1 + v,U)).*besselk(v,W) - U.*besselj(v,U).*(besselk(-1 + v,W) + besselk(1 + v,W))).*                 (2.*besselk(v,W).*(besselk(-1 + v,W) + besselk(1 + v,W)) - W.*(besselk(-1 + v,W) + besselk(1 + v,W)).^2 +                 W.*besselk(v,W).*(besselk(-2 + v,W) + 2.*besselk(v,W) + besselk(2 + v,W))) +                 U.*W.^2.*(n1.^2.*W.*(besselj(-1 + v,U) - besselj(1 + v,U)).*besselk(v,W) - n2.^2.*U.*besselj(v,U).*(besselk(-1 + v,W) + besselk(1 + v,W))).*                 (2.*besselk(v,W).*(besselk(-1 + v,W) + besselk(1 + v,W)) - W.*(besselk(-1 + v,W) + besselk(1 + v,W)).^2 +                 W.*besselk(v,W).*(besselk(-2 + v,W) + 2.*besselk(v,W) + besselk(2 + v,W)))))./(8..*n1.^2.*sqrt(-n2.^2 + neff.^2).*U.^2.*W.^5.*besselj(v,U).*besselk(v,W).^3);
        end

        function res = eigenEq_TE0m(~, ~, U, W, ~, ~, ~)
            res = (besselj(1, U)./(U.*besselj(0, U)) + besselk(1, W)./(W.*besselk(0, W)))./besselj(0, U);
        end

        function res = dev_eigenEq_TE0m(~, ~, U, W, n1, n2, k0, neff, r)
            res = -(k0.*neff.*r.*(U.*besselj(0,U).^2 + 2.*U.*besselj(1,U).^2 - besselj(0,U).*(2.*besselj(1,U) + U.*besselj(2,U))))./(2.*sqrt((n1 - neff).*(n1 + neff)).*U.^2.*besselj(0,U).^2) - (k0.*neff.*r.*(W.*besselk(0,W).^2 - 2.*W.*besselk(1,W).^2 + besselk(0,W).*(2.*besselk(1,W) + W.*besselk(2,W))))./(2.*sqrt(-n2.^2 + neff.^2).*W.^2.*besselk(0,W).^2);
            res = -res;
        end

        function res = eigenEq_TM0m(~, ~, U, W, n1, n2, ~)
            res = (n1.^2.*besselj(1, U)./(U.*besselj(0, U)) + n2.^2.*besselk(1, W)./(W.*besselk(0, W)))./besselj(0, U);
        end

        function res = dev_eigenEq_TM0m(~, ~, U, W, n1, n2, k0, neff, r)
            res = -(k0.*neff.*r.*(U.*besselj(0,U).^2 + 2.*U.*besselj(1,U).^2 - besselj(0,U).*(2.*besselj(1,U) + U.*besselj(2,U))))./(2..*sqrt((n1 - neff).*(n1 + neff)).*U.^2.*besselj(0,U).^2) - (k0.*n2.^2.*neff.*r.*(W.*besselk(0,W).^2 - 2.*W.*besselk(1,W).^2 + besselk(0,W).*(2.*besselk(1,W) + W.*besselk(2,W))))./(2..*n1.^2.*sqrt(-n2.^2 + neff.^2).*W.^2.*besselk(0,W).^2);
            res = -res;
        end

        function [V, U, W] = VUW_num(~, paras, neff)
            [k0, r, n1, n2] = deal(paras{:});
            V = k0.*r.*sqrt(n1.^2-n2.^2);
            U = real(k0.*r.*sqrt(n1.^2-neff.^2));
            W = real(k0.*r.*sqrt(neff.^2-n2.^2));
        end

        function res = calcFiberMode(obj, paras, neff, mode, varargin)
            p = inputParser;
            addParameter(p, 'v', nan);
            parse(p, varargin{:});
            
            [k0, ~, n1, n2] = deal(paras{:});
            [V, U, W] = obj.VUW_num(paras, neff);
            if isnan(p.Results.v)
                res = mode(V, U, W, n1, n2, k0);
            else
                res = mode(V, U, W, n1, n2, k0, p.Results.v);
            end
        end

        function res = calcDevFiberMode(obj, paras, neff, mode, varargin)
            p = inputParser;
            addParameter(p, 'v', nan);
            parse(p, varargin{:});
            
            [k0, r, n1, n2] = deal(paras{:});
            [V, U, W] = obj.VUW_num(paras, neff);
            if isnan(p.Results.v)
                res = mode(V, U, W, n1, n2, k0, neff, r);
            else
                res = mode(V, U, W, n1, n2, k0, neff, r, p.Results.v);
            end
        end

        function res = fieldDefine(~)
            epsilon0 = 1/(36*pi)*10^(-9);
            mu0 = 4*pi*10^(-7);
            pre1 = (epsilon0/mu0)^(1/2);

            syms k0;
            syms V U W beta_;
            syms r n1 n2;
            syms R theta;
            syms v;
            delta = 1/2 * (1 - n2^2/n1^2);
            b1 =  1/(2*U)*(besselj(v-1, U)/besselj(v, U) - besselj(v+1, U)/besselj(v, U));
            b2 = -1/(2*W)*(besselk(v-1, W)/besselk(v, W) + besselk(v+1, W)/besselk(v, W));
            F1 = (U*W/V)^2 * (b1 + (1-2*delta)*b2)/v;
            F2 = (V/(U*W))^2 * (v / (b1+b2));
            a1 = (F2-1) / 2;
            a3 = (F1-1) / 2;
            a5 = (F1-1+2*delta) / 2;
            a2 = (F2+1) / 2;
            a4 = (F1+1) / 2;
            a6 = (F1+1-2*delta) / 2;
            f = piecewise((v - 2*floor(v/2) == 1), sin(v*theta), (v - 2*floor(v/2) == 0),  cos(v*theta));
            g = piecewise((v - 2*floor(v/2) == 1), cos(v*theta), (v - 2*floor(v/2) == 0), -sin(v*theta));
            g2 = cos(2*v*theta);
            pre2 = k0*n1^2/beta_; %#ok<*CPROP>
            res = struct(...
                'HE', struct(...
                    'E', struct(...
                        'r', struct(...
                            'in',  -       (a1*besselj(v-1, U*R) + a2*besselj(v+1, U*R)) / besselj(v, U) * f,...
                            'out', - U/W * (a1*besselk(v-1, W*R) - a2*besselk(v+1, W*R)) / besselk(v, W) * f...
                        ),...
                        'phi', struct(...
                            'in',  -       (a1*besselj(v-1, U*R) - a2*besselj(v+1, U*R)) / besselj(v, U) * g,...
                            'out', - U/W * (a1*besselk(v-1, W*R) + a2*besselk(v+1, W*R)) / besselk(v, W) * g...
                        ),...
                        'z', struct(...
                            'in',  - 1j*U/(r*beta_) * besselj(v, U*R)/besselj(v, U) * f,...
                            'out', - 1j*U/(r*beta_) * besselk(v, W*R)/besselk(v, W) * f...
                        )...
                    ), ...
                    'H', struct(...
                        'r', struct(...
                            'in',  pre1 * pre2 *       (a3*besselj(v-1, U*R) - a4*besselj(v+1, U*R)) / besselj(v, U) * g,...
                            'out', pre1 * pre2 * U/W * (a5*besselk(v-1, W*R) + a6*besselk(v+1, W*R)) / besselk(v, W) * g...
                        ),...
                        'phi', struct(...
                            'in',  - pre1 * pre2 *       (a3*besselj(v-1, U*R) + a4*besselj(v+1, U*R)) / besselj(v, U) * f,...
                            'out', - pre1 * pre2 * U/W * (a5*besselk(v-1, W*R) - a6*besselk(v+1, W*R)) / besselk(v, W) * f...
                        ),...
                        'z', struct(...
                            'in',  - 1j * pre1 * U*F2/(k0*r) * besselj(v, U*R)/besselj(v, U) * g,...
                            'out', - 1j * pre1 * U*F2/(k0*r) * besselk(v, W*R)/besselk(v, W) * g...
                        )...
                    ),...
                    'Sz', struct(...
                        'in',  1/2*pre1*pre2 / besselj(v, U)^2 *            (a1*a3*besselj(v-1, U*R)^2 + a2*a4*besselj(v+1, U*R)^2 +         (1-F1*F2)/2*besselj(v-1, U*R)*besselj(v+1, U*R)*g2),...
                        'out', 1/2*pre1*pre2 / besselk(v, W)^2 * (U/W)^2 * (a1*a5*besselk(v-1, W*R)^2 + a2*a6*besselk(v+1, W*R)^2 - (1-2*delta-F1*F2)/2*besselk(v-1, W*R)*besselk(v+1, W*R)*g2)...
                    ),...
                    'Sz_int_phi_pre', struct(...
                        'in',  2*pi*r*R * 1/2*pre1*pre2 / besselj(v, U)^2 *            (a1*a3*besselj(v-1, U*R)^2 + a2*a4*besselj(v+1, U*R)^2),...
                        'out', 2*pi*r*R * 1/2*pre1*pre2 / besselk(v, W)^2 * (U/W)^2 * (a1*a5*besselk(v-1, W*R)^2 + a2*a6*besselk(v+1, W*R)^2)...
                    ),...
                    'Szintpre', struct(...
                        'in',   pi*r^2/2*pre1*pre2/besselj(v, U)^2 *            (a1*a3*                                                               R^2*(besselj(v-1, U*R)^2-besselj(v, U*R)*besselj(v-2, U*R))  + a2*a4*                                                               R^2*(besselj(v+1, U*R)^2-besselj(v, U*R)*besselj(v+2, U*R))),...
                        'out', -pi*r^2/2*pre1*pre2/besselk(v, W)^2 * (U/W)^2 * (a1*a5*((besselk(v-1, W)^2-besselk(v, W)*besselk(v-2, W)) + R^2*(-besselk(v-1, W*R)^2+besselk(v, W*R)*besselk(v-2, W*R))) + a2*a6*((besselk(v+1, W)^2-besselk(v, W)*besselk(v+2, W)) + R^2*(-besselk(v+1, W*R)^2+besselk(v, W*R)*besselk(v+2, W*R))))...
                    ),...
                    'N', struct(...
                        'in',   pi*r^2/2*pre1*pre2/besselj(v, U)^2 *            (a1*a3*(besselj(v-1, U)^2-besselj(v, U)*besselj(v-2, U)) + a2*a4*(besselj(v+1, U)^2-besselj(v, U)*besselj(v+2, U))),...
                        'out', -pi*r^2/2*pre1*pre2/besselk(v, W)^2 * (U/W)^2 * (a1*a5*(besselk(v-1, W)^2-besselk(v, W)*besselk(v-2, W)) + a2*a6*(besselk(v+1, W)^2-besselk(v, W)*besselk(v+2, W)))...
                    )...
                ),...
                'EH', struct(...
                    'E', struct(...
                        'r', struct(...
                            'in',  -       (a1*besselj(v-1, U*R) + a2*besselj(v+1, U*R)) / besselj(v, U) * f,...
                            'out', - U/W * (a1*besselk(v-1, W*R) - a2*besselk(v+1, W*R)) / besselk(v, W) * f...
                        ),...
                        'phi', struct(...
                            'in',  -       (a1*besselj(v-1, U*R) - a2*besselj(v+1, U*R)) / besselj(v, U) * g,...
                            'out', - U/W * (a1*besselk(v-1, W*R) + a2*besselk(v+1, W*R)) / besselk(v, W) * g...
                        ),...
                        'z', struct(...
                            'in',  - 1j*U/(r*beta_) * besselj(v, U*R)/besselj(v, U) * f,...
                            'out', - 1j*U/(r*beta_) * besselk(v, W*R)/besselk(v, W) * f...
                        )...
                    ), ...
                    'H', struct(...
                        'r', struct(...
                            'in',  pre1 * pre2 *       (a3*besselj(v-1, U*R) - a4*besselj(v+1, U*R)) / besselj(v, U) * g,...
                            'out', pre1 * pre2 * U/W * (a5*besselk(v-1, W*R) + a6*besselk(v+1, W*R)) / besselk(v, W) * g...
                        ),...
                        'phi', struct(...
                            'in',  - pre1 * pre2 *       (a3*besselj(v-1, U*R) + a4*besselj(v+1, U*R)) / besselj(v, U) * f,...
                            'out', - pre1 * pre2 * U/W * (a5*besselk(v-1, W*R) - a6*besselk(v+1, W*R)) / besselk(v, W) * f...
                        ),...
                        'z', struct(...
                            'in',  - 1j * pre1 * U*F2/(k0*r) * besselj(v, U*R)/besselj(v, U) * g,...
                            'out', - 1j * pre1 * U*F2/(k0*r) * besselk(v, W*R)/besselk(v, W) * g...
                        )...
                    ),...
                    'Sz', struct(...
                        'in',  1/2*pre1*pre2 / besselj(v, U)^2 *            (a1*a3*besselj(v-1, U*R)^2 + a2*a4*besselj(v+1, U*R)^2 -         (1-F1*F2)/2*besselj(v-1, U*R)*besselj(v+1, U*R)*g2),...
                        'out', 1/2*pre1*pre2 / besselk(v, W)^2 * (U/W)^2 * (a1*a5*besselk(v-1, W*R)^2 + a2*a6*besselk(v+1, W*R)^2 + (1-2*delta-F1*F2)/2*besselk(v-1, W*R)*besselk(v+1, W*R)*g2)...
                    ),...
                    'Sz_int_phi_pre', struct(...
                        'in',  2*pi*r*R * 1/2*pre1*pre2 / besselj(v, U)^2 *            (a1*a3*besselj(v-1, U*R)^2 + a2*a4*besselj(v+1, U*R)^2),...
                        'out', 2*pi*r*R * 1/2*pre1*pre2 / besselk(v, W)^2 * (U/W)^2 * (a1*a5*besselk(v-1, W*R)^2 + a2*a6*besselk(v+1, W*R)^2)...
                    ),...
                    'Szintpre', struct(...
                        'in',   pi*r^2/2*pre1*pre2/besselj(v, U)^2 *            (a1*a3*                                                               R^2*(besselj(v-1, U*R)^2-besselj(v, U*R)*besselj(v-2, U*R))  + a2*a4*                                                               R^2*(besselj(v+1, U*R)^2-besselj(v, U*R)*besselj(v+2, U*R))),...
                        'out', -pi*r^2/2*pre1*pre2/besselk(v, W)^2 * (U/W)^2 * (a1*a5*((besselk(v-1, W)^2-besselk(v, W)*besselk(v-2, W)) + R^2*(-besselk(v-1, W*R)^2+besselk(v, W*R)*besselk(v-2, W*R))) + a2*a6*((besselk(v+1, W)^2-besselk(v, W)*besselk(v+2, W)) + R^2*(-besselk(v+1, W*R)^2+besselk(v, W*R)*besselk(v+2, W*R))))...
                    ),...
                    'N', struct(...
                        'in',   pi*r^2/2*pre1*pre2/besselj(v, U)^2 *            (a1*a3*(besselj(v-1, U)^2-besselj(v, U)*besselj(v-2, U)) + a2*a4*(besselj(v+1, U)^2-besselj(v, U)*besselj(v+2, U))),...
                        'out', -pi*r^2/2*pre1*pre2/besselk(v, W)^2 * (U/W)^2 * (a1*a5*(besselk(v-1, W)^2-besselk(v, W)*besselk(v-2, W)) + a2*a6*(besselk(v+1, W)^2-besselk(v, W)*besselk(v+2, W)))...
                    )...
                ),...
                'TE', struct(...
                    'E', struct(...
                        'r', struct(...
                            'in', 0,...
                            'out', 0 ...
                        ),...
                        'phi', struct(...
                            'in',  -besselj(1, U*R)/besselj(1, U),...
                            'out', -besselk(1, W*R)/besselk(1, W)...
                        ),...
                        'z', struct(...
                            'in', 0,...
                            'out', 0 ...
                        )...
                    ),...
                    'H', struct(...
                        'r', struct(...
                            'in',  pre1 * beta_/k0 * besselj(1, U*R)/besselj(1, U),...
                            'out', pre1 * beta_/k0 * besselk(1, W*R)/besselk(1, W)...
                        ),...
                        'phi', struct(...
                            'in', 0,...
                            'out', 0 ...
                        ),...
                        'z', struct(...
                            'in',   1j * pre1 * U/(k0*r) * besselj(0, U*R)/besselj(1, U),...
                            'out', -1j * pre1 * W/(k0*r) * besselk(0, W*R)/besselk(1, W)...
                        )...
                    ),...
                    'Sz', struct(...
                        'in',  1/2*pre1*beta_/k0*besselj(1, U*R)^2/besselj(1, U)^2,...
                        'out', 1/2*pre1*beta_/k0*besselk(1, W*R)^2/besselk(1, W)^2 ...
                    ),...
                    'Sz_int_phi_pre', struct(...
                        'in',  2*pi*r*R * 1/2*pre1*beta_/k0*besselj(1, U*R)^2/besselj(1, U)^2,...
                        'out', 2*pi*r*R * 1/2*pre1*beta_/k0*besselk(1, W*R)^2/besselk(1, W)^2 ...
                    ),...
                    'Szintpre', struct(...
                        'in',   2*pi*r^2 /2*pre1*beta_/k0 * 1/2*(                                                                              R^2*(besselj(1, U*R)^2-besselj(0, U*R)*besselj(2, U*R))/besselj(1, U)^2),...
                        'out', -2*pi*r^2 /2*pre1*beta_/k0 * 1/2*((besselk(1, W)^2-besselk(0, W)*besselk(2, W))/besselk(1, W)^2 - R^2*(-besselk(1, W*R)^2+besselk(0, W*R)*besselk(2, W*R))/besselk(1, W)^2)...
                    ),...
                    'N', struct(...
                        'in',   2*pi*r^2/2*pre1*beta_/k0 * 1/2*(besselj(1, U)^2-besselj(0, U)*besselj(2, U))/besselj(1, U)^2,...
                        'out', -2*pi*r^2/2*pre1*beta_/k0 * 1/2*(besselk(1, W)^2-besselk(0, W)*besselk(2, W))/besselk(1, W)^2 ...
                    )...
                ),...
                'TM', struct(...
                    'E', struct(...
                        'r', struct(...
                            'in',               besselj(1, U*R)/besselj(1, U),...
                            'out', (n1/n2)^2 * besselk(1, W*R)/besselk(1, W)...
                        ),...
                        'phi', struct(...
                            'in', 0,...
                            'out', 0 ...
                        ),...
                        'z', struct(...
                            'in',   1j*U/(beta_*r) *              besselj(0, U*R)/besselj(1, U),...
                            'out', -1j*W/(beta_*r) * (n1/n2)^2 * besselk(0, W*R)/besselk(1, W)...
                        )...
                    ),...
                    'H', struct(...
                        'r', struct(...
                            'in', 0,...
                            'out', 0 ...
                        ),...
                        'phi', struct(...
                            'in',  pre1 * pre2 * besselj(1, U*R)/besselj(1, U),...
                            'out', pre1 * pre2 * besselk(1, W*R)/besselk(1, W)...
                        ),...
                        'z', struct(...
                            'in', 0,...
                            'out', 0 ...
                        )...
                    ),...
                    'Sz', struct(...
                        'in',  1/2*pre1*pre2*            besselj(1, U*R)^2/besselj(1, U)^2,...
                        'out', 1/2*pre1*pre2/(1-2*delta)*besselk(1, W*R)^2/besselk(1, W)^2 ...
                    ),...
                    'Sz_int_phi_pre', struct(...
                        'in',  2*pi*r*R * 1/2*pre1*pre2*            besselj(1, U*R)^2/besselj(1, U)^2,...
                        'out', 2*pi*r*R * 1/2*pre1*pre2/(1-2*delta)*besselk(1, W*R)^2/besselk(1, W)^2 ...
                    ),...
                    'Szintpre', struct(...
                        'in',   2*pi*r^2 /2*pre1*pre2*              1/2*(                                                                              R^2*(besselj(1, U*R)^2-besselj(0, U*R)*besselj(2, U*R))/besselj(1, U)^2),...
                        'out', -2*pi*r^2 /2*pre1*pre2/(1-2*delta) * 1/2*((besselk(1, W)^2-besselk(0, W)*besselk(2, W))/besselk(1, W)^2 - R^2*(-besselk(1, W*R)^2+besselk(0, W*R)*besselk(2, W*R))/besselk(1, W)^2)...
                    ),...
                    'N', struct(...
                        'in',   2*pi*r^2 /2*pre1*pre2*              1/2*(besselj(1, U)^2-besselj(0, U)*besselj(2, U))/besselj(1, U)^2,...
                        'out', -2*pi*r^2 /2*pre1*pre2/(1-2*delta) * 1/2*(besselk(1, W)^2-besselk(0, W)*besselk(2, W))/besselk(1, W)^2 ...
                    )...
                )...
            );

            for mode0 = fieldnames(res)'
                mode = mode0{1};
                res.(mode).Sztotal = res.(mode).N.in + res.(mode).N.out;
                res.(mode).Szint = struct(...
                    'in', res.(mode).Szintpre.in/res.(mode).Sztotal,...
                    'out', (res.(mode).N.in + res.(mode).Szintpre.out)/res.(mode).Sztotal...
                );
                res.(mode).Sz_int_phi = struct(...
                    'in', res.(mode).Sz_int_phi_pre.in/res.(mode).Sztotal,...
                    'out', res.(mode).Sz_int_phi_pre.out/res.(mode).Sztotal...
                );
                for EH = ['E', 'H']
                    res.(mode).(EH).x = struct(...
                        'in',  res.(mode).(EH).r.in*cos(theta) - res.(mode).(EH).phi.in*sin(theta),...
                        'out', res.(mode).(EH).r.out*cos(theta) - res.(mode).(EH).phi.out*sin(theta)...
                    );
                    res.(mode).(EH).y = struct(...
                        'in',  res.(mode).(EH).r.in*sin(theta) + res.(mode).(EH).phi.in*cos(theta),...
                        'out', res.(mode).(EH).r.out*sin(theta) + res.(mode).(EH).phi.out*cos(theta)...
                    );
                end
                normField = @(f, io) f.x.(io)*conj(f.x.(io)) + f.y.(io)*conj(f.y.(io)) + f.z.(io)*conj(f.z.(io));
                for EH = ['E', 'H']
                    res.(mode).(EH).norm = struct(...
                        'in', normField(res.(mode).(EH), 'in'),...
                        'out', normField(res.(mode).(EH), 'out')...
                    );
                end
            end
        end
        
        function res = lamd(~, x)
            if (isstruct(x) && isequal(x.in, 0)) || (~isstruct(x) && isequal(x, 0))
                res = @(V, U, W, R, theta, n1, n2, beta_, r, k0, v) zeros(size(R));
            else
                syms V U W R theta n1 n2 beta_ r k0 v;
                paras = [V, U, W, R, theta, n1, n2, beta_, r, k0, v]; %#ok<*CPROPLC>
                if isstruct(x)
                    xin = matlabFunction(x.in, 'Vars', paras, 'File', tempname('tempFunctions\'));
                    xout = matlabFunction(x.out, 'Vars', paras, 'File', tempname('tempFunctions\'));
                    res = @(V, U, W, R, theta, n1, n2, beta_, r, k0, v) (R<=1).*xin(V, U, W, R, theta, n1, n2, beta_, r, k0, v) + (R>1).*xout(V, U, W, R, theta, n1, n2, beta_, r, k0, v);
                else
                    res = matlabFunction(x, 'Vars', paras, 'File', tempname('tempFunctions\'));
                end
            end
        end

        function res0 = mergefunc_(obj, res, func, key, field, varargin)
            p = inputParser;
            addParameter(p, 'v', nan);
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;

            foofun = @foo;
            function foo_ = foo(indeps, key, delam)
                neff = res.(key).neff(indeps{:});
                temp = obj.basicParaCal(indeps, 'usage', 'foo', 'output', 'simple');
                [lam, k0, r, n1, n2] = deal(temp{:});
                lam = lam + delam;
                if ~isnumeric(obj.n1lam)
                    n1 = obj.n1lam(lam);
                end
                if ~isnumeric(obj.n2lam)
                    n2 = obj.n2lam(lam);
                end
                [V, U, W] = obj.VUW_num({k0, r, n1, n2}, neff);
                beta_ = k0.*neff;
                if ismember(field, {'Szint', 'Sz_int_phi'})
                    foo_ = @(rmesh) func(V, U, W, rmesh./r, 0, n1, n2, beta_, r, k0, p.Results.v);
                else
                    foo_ = @(rmesh, thetamesh) func(V, U, W, rmesh./r, thetamesh, n1, n2, beta_, r, k0, p.Results.v);
                end
            end
            if obj.ndim == 0
                res0 = @() foofun({}, key, delam);
            else
                res0 = @(x) foofun({x}, key, delam);
            end
        end

        function res0 = eta(obj, res, key, indeps, varargin)
            p = inputParser;
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;

            temp = res.(key).Szint(indeps{:});
            temp2 = obj.basicParaCal(indeps, 'usage', 'eta', 'output', 'simple', 'delam', delam);
            [~, ~, r, ~, ~] = deal(temp2{:});
            res0 = temp(r);
        end

        function res0 = Reff(obj, res, key, indeps, varargin)
            p = inputParser;
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;

            Szintfun = res.(key).Szint(indeps{:});
            Szintfunprime = res.(key).Sz_int_phi(indeps{:});
            tempfun = @(r) deal(Szintfun(r) - (1-1/exp(1)^2), Szintfunprime(r));

            temp = obj.basicParaCal(indeps, 'usage', 'Reff', 'output', 'simple', 'delam', delam);
            [~, ~, r, ~, ~] = deal(temp{:});

            options = optimoptions('fsolve','Display','off','SpecifyObjectiveGradient',true, 'FunctionTolerance', 1e-10, 'OptimalityTolerance', 1e-200);
            res0 = nan;
            for ii=[0.9, 10, 20, 50, 1e2, 200, 500, 1e3, 2e3, 5e3, 1e4]
                [res00, fval] = fsolve(tempfun, ii*r, options);
                if fval < 1e-6 && res00 > 0
                    res0 = res00;
                    break
                end
            end
        end

        function res0 = vg(obj, res, key, indeps, varargin)
            p = inputParser;
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;

            neff = res.(key).neff(indeps{:});
            temp = obj.basicParaCal(indeps, 'output', 'simple', 'delam', delam, 'usage', 'vg');
            [~, ~, ~, n1, n2] = deal(temp{:});
            delta = 1/2 .* (1 - n2.^2./n1.^2);
            c = 3e8;
            eta = res.(key).eta(indeps{:});
            res0 = c./n1.^2.*neff.*1./(1-2.*delta.*(1-eta));
        end

        function foo = deco(obj, func)
            function res = foo_(points)
                res = [];
                for i = 1:length(points)
                    res(i) = func(points(i));
                end
            end
            if obj.ndim == 0
                foo = func;
            else
                foo = @foo_;   
            end
        end
            
        function res = fieldCalc(obj, res, varargin)
            p = inputParser;
            addParameter(p, 'delam', 0);
            parse(p, varargin{:});
            delam = p.Results.delam;
            
            for key0 = fieldnames(res)'
                key = key0{1};
                mode = key(1:2);
                v = nan;
                if ismember(mode, {'HE', 'EH'})
                    v = str2double(key(3));
                end
                for field0 = {'E', 'H', 'Sz', 'Szint', 'Sz_int_phi'}
                    field = field0{1};
                    mergefunc = @(fun) obj.mergefunc_(res, fun, key, field, 'v', v, 'delam', delam);
                    comp = @(fun) mergefunc(obj.lamd(fun));
                    if ismember(field, {'E', 'H'})
                        res.(key).(field) = structfun(comp, obj.fields.(mode).(field), 'UniformOutput', false);
                    elseif ismember(field, {'Sz', 'Szint', 'Sz_int_phi'})
                        res.(key).(field) = comp(obj.fields.(mode).(field));
                    end
                end
                
                if obj.ndim == 0
                    res.(key).eta = obj.deco(@() obj.eta(res, key, {}, 'delam', delam));
                    res.(key).Reff = obj.deco(@() obj.Reff(res, key, {}, 'delam', delam));
                    res.(key).vg = obj.deco(@() obj.vg(res, key, {}, 'delam', delam));
                else
                    res.(key).eta = obj.deco(@(indep) obj.eta(res, key, {indep}, 'delam', delam));
                    res.(key).Reff = obj.deco(@(indep) obj.Reff(res, key, {indep}, 'delam', delam));
                    res.(key).vg = obj.deco(@(indep) obj.vg(res, key, {indep}, 'delam', delam));
                end
            end  
        end

        function obj = DwCalc(obj)
            dedis = {};
            delams = [1e-3, 2e-3, 5e-3, 1e-2]*1e-9;
            for delam = delams
                dedis0 = {};
                for sign = [1, -1]
                    temp = obj.basicParaCal(obj.indeps_vals, 'neff', obj.neffmlist_, 'delam', sign*delam, 'output', 'simple');
                    [~, k0, r, n1, n2] = deal(temp{:}); %#ok<*PROP>
                    dis = obj.dispersiveCurveCalc({k0, r, n1, n2}, 'delam', sign*delam);
                    dis = obj.fieldCalc(dis, 'delam', sign*delam);
                    dedis0{end+1} = dis; %#ok<*AGROW>
                end
                dedis{end+1} = dedis0;
            end
            for key0 = fieldnames(obj.res)'
                key = key0{1};
                indeps_Dw = [];
                if obj.ndim == 0
                    vg_0 = obj.res.(key).vg();
                    vg_pre2 = cellfun(@(x) cellfun(@(y) y.(key).vg(), x), dedis, 'UniformOutput', false);
                    vg_pre = [];
                    for j = 1:length(vg_pre2)
                        for k = 1:length(vg_pre2{1})
                            vg_pre(j, k) = vg_pre2{j}(k);
                        end
                    end
                    Dw_pre = ((1./vg_pre-1./vg_0)'./[delams; -delams]);
                    Dw_pre = mean(Dw_pre);
                    Dw = interp1(delams, Dw_pre, 0, 'spline', 'extrap');
                    obj.res.(key).Dw = @() Dw;
                else
                    for i = 1:length(obj.indeps_vals)
                        indep = obj.indeps_vals(i);
                        vg_0 = obj.res.(key).vg(indep);
                        if isnan(vg_0)
                            continue
                        end
                        vg_pre2 = cellfun(@(x) cellfun(@(y) y.(key).vg(indep), x), dedis, 'UniformOutput', false);
                        vg_pre = [];
                        for j = 1:length(vg_pre2)
                            for k = 1:length(vg_pre2{1})
                                vg_pre(j, k) = vg_pre2{j}(k);
                            end
                        end
                        Dw_pre = ((1./vg_pre-1./vg_0)'./[delams; -delams]);
                        Dw_pre = mean(Dw_pre);
                        if any(isnan(Dw_pre))
                            continue
                        end
                        Dw = interp1(delams, Dw_pre, 0, 'spline', 'extrap');
                        if ~isnan(Dw)
                            indeps_Dw(i, :) = [indep, Dw];
                        end
                    end
                    indeps_Dw = indeps_Dw';
                    obj.res.(key).Dw = @(x) interp1(indeps_Dw(1, :), indeps_Dw(2, :), x, 'spline', nan);
                end
            end
        end
    end
end
