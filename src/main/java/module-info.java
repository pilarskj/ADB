open module adb {
    requires beast.pkgmgmt;
    requires beast.base;
    requires org.apache.commons.statistics.distribution;
    requires commons.math3;

    exports adb;

    provides beast.base.core.BEASTInterface with
        adb.GammaBranchingModel;
}
