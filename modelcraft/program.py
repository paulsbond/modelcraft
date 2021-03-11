class Program:
    def __init__(self, executable: str):
        self.executable = executable
        self.args = []
        self.stdin = []
        self.environ = {}

    def script(self) -> str:
        script = "#!/usr/bin/env bash\n\n"
        if self.environ:
            for variable, value in self.environ.items():
                script += f"export {variable}={value}\n"
            script += "\n"
        script += self.executable
        script += f" {' '.join(self.args)} \\\n> stdout.txt 2> stderr.txt"
        if self.stdin:
            script += " << EOF\n"
            for line in self.stdin:
                script += f"{line}\n"
            script += "EOF\n"
        else:
            script += "\n"
        return script
