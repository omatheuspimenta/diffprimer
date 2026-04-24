from typing import Annotated

from typer import Option, Typer

app = Typer()

@app.command()
def main(
    version: Annotated[bool, Option("--version", is_eager=True)] = False,
    reference_file: str = Option(...),
):
    print("Running with", reference_file)

if __name__ == "__main__":
    app()
