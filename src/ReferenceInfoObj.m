classdef ReferenceInfoObj

    methods (Static = true)
        function [javaInstance_ret, build_ret] = getJavaInstance(refdir)
            persistent javaInstance build
            if nargin >0
                javaInstance = org.broadinstitute.cga.tools.seq.ReferenceInfo(refdir);
                build = javaInstance.getBuild();
            end
            javaInstance_ret = javaInstance;
            build_ret = build;
        end
        function check_build(javaInstance, build)
            if ~strcmp(build, javaInstance.getBuild())
                error (['input build ' build ' not initialized: ' char(javaInstance.getBuild()) ]);
                %error ('only one build supported now');
            end
        end
        function initialized = isInitialized()
           initialized = true;
           try
               java_instance = ReferenceInfoObj.getJavaInstance();
               if isempty(java_instance)
                   initialized = false;
               end
           catch ME
               initialized = false;
           end
        end
        
        function init (refdir)
            ReferenceInfoObj.getJavaInstance(refdir);
            % no return value
        end
        
        function chrom = getChrom(num, build)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            ReferenceInfoObj.check_build(javaInstance, build);
            chrom = cell(javaInstance.getChrom(num));
            
        end
        
        function num = getNum(chrom, build)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            ReferenceInfoObj.check_build(javaInstance, build);
            num = javaInstance.getNum(chrom);
            
        end
        
        function len = getLength(num, build)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            ReferenceInfoObj.check_build(javaInstance, build);
            if ~isnumeric(num)
                num=javaInstance.getNum(num);
            end
            len = javaInstance.getLength(num);
            
        end

        
        function maxnum = getMaxNum(build)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            ReferenceInfoObj.check_build(javaInstance, build);
            maxnum = javaInstance.getMaxNum();            
        end
        
        function use = getUse(num,build)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            ReferenceInfoObj.check_build(javaInstance, build);
            if ~isnumeric(num)
                num=javaInstance.getNum(num);
            end
            use = cell(javaInstance.getUse(num));
        end

        
        function MaleCN = getMaleCN(num, build)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            ReferenceInfoObj.check_build(javaInstance, build);
            if ~isnumeric(num)
                num=javaInstance.getNum(num);
            end
            MaleCN = javaInstance.getMaleCN(num);
            
        end

        
        function FemaleCN = getFemaleCN(num, build)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            ReferenceInfoObj.check_build(javaInstance, build);
            if ~isnumeric(num)
                num=javaInstance.getNum(num);
            end                
            FemaleCN = javaInstance.getFemaleCN(num);            
        end 

        
        function build = getBuild()
            javaInstance = ReferenceInfoObj.getJavaInstance();
            build = char(javaInstance.getBuild());            
        end
        
        function refFile = getReferenceFile(num)
            javaInstance = ReferenceInfoObj.getJavaInstance();
            refFile = char(javaInstance.getReferenceFile(num));            
        end
    end    
end
