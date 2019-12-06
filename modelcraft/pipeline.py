#!/usr/bin/python3

import modelcraft.args as args
import modelcraft.gemmineer as gemmineer
import modelcraft.jobs as jobs
import modelcraft.report as report
import shutil


class Pipeline:
    def __init__(self, arguments):
        print("# ModelCraft")
        print("\nPlease cite [paper to be published]")
        args.parse(arguments)
        self.min_rwork = 1
        self.min_rfree = 1
        self.min_fragments_built = 999
        self.max_longest_fragment = 0
        self.max_residues_built = 0
        self.max_residues_sequenced = 0
        self.cycles_without_improvement = 0
        self.run()

    def run(self):
        for self.cycle in range(1, args.cycles + 1):
            print("\n## Cycle %d\n" % self.cycle)
            mlhl = args.unbiased and self.min_rwork > 0.35
            if self.cycle == 1:
                missingPhases = args.colin_hl is None and args.colin_phifom is None
                if missingPhases and args.mr_model is not None:
                    refmac = jobs.Refmac(self.cycle, args.mr_model, mlhl)
                    buccaneer = jobs.Buccaneer(self.cycle, refmac.hklout, args.xyzin, cycles=3)
                else:
                    buccaneer = jobs.Buccaneer(self.cycle, args.mtzin, args.xyzin, cycles=3)
            else:
                coot = jobs.Prune(self.cycle, refmac.xyzout, refmac.hklout)
                refmac = jobs.Refmac(self.cycle, coot.xyzout, mlhl)
                buccaneer = jobs.Buccaneer(self.cycle, refmac.hklout, refmac.xyzout)
            refmac = jobs.Refmac(self.cycle, buccaneer.xyzout, mlhl)
            coot = jobs.Prune(self.cycle, refmac.xyzout, refmac.hklout, chains_only=True)
            if args.add_waters and self.min_rwork < 0.4:
                coot = jobs.FindWaters(self.cycle, coot.xyzout, refmac.hklout)
            refmac = jobs.Refmac(self.cycle, coot.xyzout, mlhl)
            self.process_cycle_output(refmac)
            if self.cycle > 1:
                jobs.remove_job_directories(self.cycle - 1)
            if args.auto_stop and self.cycles_without_improvement == 4:
                break
        jobs.remove_job_directories(self.cycle)

    def improved(self, cycle):
        required_improvement = 0.02
        improvement = (self.min_rwork - cycle["rwork"]) / self.min_rwork
        if improvement > required_improvement:
            return True
        improvement = (cycle["residues_built"] - self.max_residues_built) / float(self.max_residues_built)
        if improvement > required_improvement:
            return True
        improvement = (cycle["residues_sequenced"] - self.max_residues_sequenced) / float(self.max_residues_sequenced)
        if improvement > required_improvement:
            return True
        improvement = (self.min_fragments_built - cycle["fragments_built"]) / float(self.min_fragments_built)
        if improvement > required_improvement:
            return True
        improvement = (cycle["longest_fragment"] - self.max_longest_fragment) / float(self.max_longest_fragment)
        if improvement > required_improvement:
            return True
        return False

    def process_cycle_output(self, refmac):
        cycle = {
            "cycle": self.cycle,
            "rwork": refmac.final_rwork,
            "rfree": refmac.final_rfree,
        }
        cycle.update(gemmineer.model_stats(refmac.xyzout))
        print("\nResidues built: %d" % cycle["residues_built"])
        print("Residues sequenced: %d" % cycle["residues_sequenced"])
        print("R-work: %.3f" % cycle["rwork"])
        print("R-free: %.3f" % cycle["rfree"])

        if self.improved(cycle):
            self.cycles_without_improvement = 0
        else:
            self.cycles_without_improvement += 1
            print("\nNo significant improvement for %d cycle(s)" % self.cycles_without_improvement)

        if cycle["rwork"] < self.min_rwork:
            self.min_rwork = cycle["rwork"]
        if cycle["rfree"] < self.min_rfree:
            self.min_rfree = cycle["rfree"]
            print("Copying files to output as R-free improved")
            shutil.copyfile(str(refmac.xyzout), "xyzout.pdb")
            shutil.copyfile(str(refmac.hklout), "hklout.mtz")
        if cycle["residues_built"] > self.max_residues_built:
            self.max_residues_built = cycle["residues_built"]
        if cycle["residues_sequenced"] > self.max_residues_sequenced:
            self.max_residues_sequenced = cycle["residues_sequenced"]
        if cycle["fragments_built"] < self.min_fragments_built:
            self.min_fragments_built = cycle["fragments_built"]
        if cycle["longest_fragment"] > self.max_longest_fragment:
            self.max_longest_fragment = cycle["longest_fragment"]

        report.add_cycle(cycle)
