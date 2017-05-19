classdef Pipe     
    properties (Access = private)
        arg
    end

      methods
          function this = Pipe(arg)
              this.arg = arg;
          end

          function result = subsref(this, s)
              if mod(numel(s), 2) ~= 0 || any(~strcmp({s(1:2:end).type}, '.')) || any(~strcmp({s(2:2:end).type}, '()'))
                  error('invalid syntax');
              end
              result = this.arg;
              argpos = 1;
              for idx = 1:2:numel(s)
                  if strcmp(s(idx).subs, 'as')
                      argpos = s(idx+1).subs{1};
                  else
                      args = [s(idx+1).subs(1:argpos-1), {result}, s(idx+1).subs(argpos:end)];
                      result = feval(s(idx).subs, args{:});
                      argpos = 1;
                  end
              end                
          end
      end    
  end
