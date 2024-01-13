classdef NuPrinter < handle
    %NUPRINTER This class is nuclear printer object defination

    properties(SetAccess=private, GetAccess=private, Hidden)
        caller
        ax_caller

        outargs
        slice_n
        fname
    end

    properties(Constant, Hidden)
        STATUS_SUCCESS       =   0
        STATUS_WRITE_ERROR   =   -1
    end
    
    methods(Access=public)
        function this = NuPrinter(volopts, outargs)
            arguments
                volopts (1,1) struct
                outargs (1,1) struct
            end

            this.outargs = outargs;
            this.slice_n = volopts.slices;
            % file name without suffix
            this.fname = outargs.folder+filesep+outargs.file;
        end

        function SetCaller(this, caller)
            this.caller = caller;
            this.ax_caller = caller.UIAxes_Viewer;
        end
        
        function [status, file] = print(this, direct)
            %PRINT This function prints objects and redirect to files
            arguments
                this
                direct (1,1) string {mustBeMember(direct, ...
                    ["normal", "reverse"])} = "normal";
            end

            switch direct
                case "normal"
                    this.caller.MoveSliceHome();
                    slices_range = 1:this.slice_n;
                case "reverse"
                    this.caller.MoveSliceEnd();
                    slices_range = this.slice_n:-1:1;
            end

            this.caller.SetSystemInfo("export results...");

            switch this.outargs.format
                case {'gif','tif'}   % using imwrite

                case "avi"           % using VideoWriter

                case "png"           % using exportgraphics
                    for pslice = slices_range
                        file = sprintf("%s_%d.png", this.fname, pslice);
                        exportgraphics(this.ax_caller, file, ...
                            "Resolution",this.outargs.resolution);
                        if direct == "normal"
                            this.caller.SetProgressBar(pslice/this.slice_n);
                            this.caller.MoveSliceNext();
                        else
                            this.caller.SetProgressBar((this.slice_n-pslice+1)/this.slice_n);
                            this.caller.MoveSlicePrev();
                        end
                    end
            end

            this.caller.SetSystemInfo("export success");
            status = this.STATUS_SUCCESS;
        end
    end

    methods(Access=private)
        function normal_print(this)

        end

        function reverse_print(this)

        end
    end
end

