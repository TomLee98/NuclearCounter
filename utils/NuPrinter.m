classdef NuPrinter < handle
    %NUPRINTER This class is nuclear printer object defination

    properties(SetAccess=private, GetAccess=private, Hidden)
        caller
        ax_caller
    end

    properties(SetAccess=private, GetAccess=private)
        outarg
        slice_n
        fname
    end

    properties(Constant, Hidden)
        STATUS_SUCCESS       =   0
        STATUS_WRITE_ERROR   =   -1
    end
    
    methods(Access=public)
        function this = NuPrinter(vopt, oarg)
            arguments
                vopt (1,12) table
                oarg  (1,1)  struct
            end

            this.outarg = oarg;
            this.slice_n = vopt.slices;
            % file name without suffix
            this.fname = string([oarg.folder, filesep])+oarg.file;
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

            this.caller.SetSystemLog("export results...");

            switch this.outarg.format
                case {'gif','tif'}   % using imwrite
                    % copy object to hidden figure for fast saving
                    pos = this.caller.GetViewerPosition();
                    fig = figure("Position",pos, ...
                                 "Visible","off");
                    if this.outarg.format == "gif"
                        file = sprintf("%s_labeled.gif", this.fname);
                        for pslice = slices_range
                            clf(fig);
                            copyobj(this.ax_caller, fig);
                            frame = getframe(gca);
                            img = frame2im(frame);
                            [A, map] = rgb2ind(img, 256);
                            if direct == "normal"
                                if pslice == 1
                                    imwrite(A,map,file,"gif","LoopCount",inf,...
                                        "DelayTime",1/this.outarg.fr);
                                else
                                    imwrite(A,map,file,"gif","WriteMode","append",...
                                        "DelayTime",1/this.outarg.fr);
                                end
                                this.caller.SetProgressBar(pslice/this.slice_n);
                                this.caller.MoveSliceNext();
                            else
                                if pslice == this.slice_n
                                    imwrite(A,map,file,"gif","LoopCount",inf,...
                                        "DelayTime",1/this.outarg.fr);
                                else
                                    imwrite(A,map,file,"gif","WriteMode","append",...
                                        "DelayTime",1/this.outarg.fr);
                                end
                                this.caller.SetProgressBar((this.slice_n-pslice+1)/this.slice_n);
                                this.caller.MoveSlicePrev();
                            end
                        end
                    else
                        file = sprintf("%s_labeled.tif", this.fname);
                        for pslice = slices_range
                            clf(fig);
                            copyobj(this.ax_caller, fig);
                            frame = getframe(gca);
                            img = frame2im(frame);
                            if direct == "normal"
                                if pslice == 1
                                    imwrite(img,file,"tif", ...
                                        "Resolution",this.outarg.resolution);
                                else
                                    imwrite(img,file,"tif","WriteMode","append", ...
                                        "Resolution",this.outarg.resolution);
                                end
                                this.caller.SetProgressBar(pslice/this.slice_n);
                                this.caller.MoveSliceNext();
                            else
                                if pslice == this.slice_n
                                    imwrite(img,file,"tif", ...
                                        "Resolution",this.outarg.resolution);
                                else
                                    imwrite(img,file,"tif","WriteMode","append", ...
                                        "Resolution",this.outarg.resolution);
                                end
                                this.caller.SetProgressBar((this.slice_n-pslice+1)/this.slice_n);
                                this.caller.MoveSlicePrev();
                            end
                        end
                    end

                    close(fig);
                case "avi"           % using VideoWriter
                    % copy object to hidden figure for fast saving
                    pos = this.caller.GetViewerPosition();
                    fig = figure("Position",pos, ...
                                 "Visible","off");
                    file = sprintf("%s_labeled.avi", this.fname);

                    vw = VideoWriter(file, "MPEG-4");
                    vw.Quality = 100;
                    vw.FrameRate = this.outarg.fr;
                    open(vw);
                    for pslice = slices_range
                        clf(fig);
                        copyobj(this.ax_caller, fig);
                        frame = getframe(gca);
                        writeVideo(vw, frame);
                        if direct == "normal"
                            this.caller.SetProgressBar(pslice/this.slice_n);
                            this.caller.MoveSliceNext();
                        else
                            this.caller.SetProgressBar((this.slice_n-pslice+1)/this.slice_n);
                            this.caller.MoveSlicePrev();
                        end
                    end
                    close(vw);
                case "png"           % using exportgraphics
                    for pslice = slices_range
                        file = sprintf("%s_labeled_%d.png", this.fname, pslice);
                        exportgraphics(this.ax_caller, file, ...
                            "Resolution",this.outarg.resolution);
                        if direct == "normal"
                            this.caller.SetProgressBar(pslice/this.slice_n);
                            this.caller.MoveSliceNext();
                        else
                            this.caller.SetProgressBar((this.slice_n-pslice+1)/this.slice_n);
                            this.caller.MoveSlicePrev();
                        end
                    end
            end

            this.caller.SetSystemLog("export success");
            status = this.STATUS_SUCCESS;
        end
    end
end

