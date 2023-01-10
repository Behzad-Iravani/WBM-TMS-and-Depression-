classdef assess_free_runing_frq

    properties
        intrinsic
        roi_name
    end
    methods
        function as_fq = assess_free_runing_frq
            as_fq.intrinsic = [0.0253906250000000
                0.0253906250000000
                0.0195312500000000
                0.0234375000000000
                0.0351562500000000
                0.0195312500000000
                0.0234375000000000
                0.0195312500000000
                0.0234375000000000
                0.0175781250000000
                0.0253906250000000
                0.0234375000000000
                0.0214843750000000
                0.0234375000000000
                0.0292968750000000
                0.0253906250000000
                0.0253906250000000
                0.0292968750000000
                0.0273437500000000
                0.0234375000000000
                0.0234375000000000
                0.0234375000000000
                0.0214843750000000
                0.0273437500000000
                0.0156250000000000
                0.0234375000000000
                0.0273437500000000
                0.0195312500000000
                0.0234375000000000
                0.0195312500000000
                0.0195312500000000
                0.0234375000000000
                0.0234375000000000
                0.0195312500000000
                0.0214843750000000
                0.0234375000000000
                0.0273437500000000
                0.0234375000000000
                0.0351562500000000
                0.0214843750000000
                0.0234375000000000
                0.0175781250000000
                0.0234375000000000
                0.0195312500000000
                0.0253906250000000
                0.0195312500000000
                0.0195312500000000
                0.0195312500000000
                0.0332031250000000
                0.0253906250000000
                0.0156250000000000
                0.0273437500000000
                0.0195312500000000
                0.0273437500000000
                0.0234375000000000
                0.0292968750000000
                0.0234375000000000
                0.0253906250000000
                0.0175781250000000
                0.0156250000000000
                0.0214843750000000
                0.0175781250000000
                0.0195312500000000
                0.0156250000000000
                0.0234375000000000
                0.0253906250000000
                0.0273437500000000
                0.0195312500000000]; % Estimated intrinsic frequency from empirical data

            roi_name = {"ctx-lh-bankssts"
                "ctx-lh-caudalanteriorcingulate"
                "ctx-lh-caudalmiddlefrontal"
                "ctx-lh-cuneus"
                "ctx-lh-entorhinal"
                "ctx-lh-fusiform"
                "ctx-lh-inferiorparietal"
                "ctx-lh-inferiortemporal"
                "ctx-lh-isthmuscingulate"
                "ctx-lh-lateraloccipital"
                "ctx-lh-lateralorbitofrontal"
                "ctx-lh-lingual"
                "ctx-lh-medialorbitofrontal"
                "ctx-lh-middletemporal"
                "ctx-lh-parahippocampal"
                "ctx-lh-paracentral"
                "ctx-lh-parsopercularis"
                "ctx-lh-parsorbitalis"
                "ctx-lh-parstriangularis"
                "ctx-lh-pericalcarine"
                "ctx-lh-postcentral"
                "ctx-lh-posteriorcingulate"
                "ctx-lh-precentral"
                "ctx-lh-precuneus"
                "ctx-lh-rostralanteriorcingulate"
                "ctx-lh-rostralmiddlefrontal"
                "ctx-lh-superiorfrontal"
                "ctx-lh-superiorparietal"
                "ctx-lh-superiortemporal"
                "ctx-lh-supramarginal"
                "ctx-lh-frontalpole"
                "ctx-lh-temporalpole"
                "ctx-lh-transversetemporal"
                "ctx-lh-insula"
                "ctx-rh-bankssts"
                "ctx-rh-caudalanteriorcingulate"
                "ctx-rh-caudalmiddlefrontal"
                "ctx-rh-cuneus"
                "ctx-rh-entorhinal"
                "ctx-rh-fusiform"
                "ctx-rh-inferiorparietal"
                "ctx-rh-inferiortemporal"
                "ctx-rh-isthmuscingulate"
                "ctx-rh-lateraloccipital"
                "ctx-rh-lateralorbitofrontal"
                "ctx-rh-lingual"
                "ctx-rh-medialorbitofrontal"
                "ctx-rh-middletemporal"
                "ctx-rh-parahippocampal"
                "ctx-rh-paracentral"
                "ctx-rh-parsopercularis"
                "ctx-rh-parsorbitalis"
                "ctx-rh-parstriangularis"
                "ctx-rh-pericalcarine"
                "ctx-rh-postcentral"
                "ctx-rh-posteriorcingulate"
                "ctx-rh-precentral"
                "ctx-rh-precuneus"
                "ctx-rh-rostralanteriorcingulate"
                "ctx-rh-rostralmiddlefrontal"
                "ctx-rh-superiorfrontal"
                "ctx-rh-superiorparietal"
                "ctx-rh-superiortemporal"
                "ctx-rh-supramarginal"
                "ctx-rh-frontalpole"
                "ctx-rh-temporalpole"
                "ctx-rh-transversetemporal"
                "ctx-rh-insula"};

        end


        function plot_intrinsic_freq(obj)
            histogram(obj.intrinsic)

        end

    end
end
