import click

from algorithms.aco import aco
from algorithms.dnala import dnala
from algorithms.genetic_algorithm import genetic_algorithm


@click.group()
@click.pass_context
@click.option('-i', '--instance', required=False, type=click.Choice(['10', '20', '30', '40', '50', '60', '70', '80', '90', '100', '125', '150', '175', '200', '300']), default='10')
@click.option('-t', '--time', required=False, type=int, default=60)
@click.option('-o', '--output', required=False, type=str, default='output.txt')
def cli(ctx, instance, time, output):
    # ensure that ctx.obj exists and is a dict (in case `cli()` is called
    # by means other than the `if` block below)
    ctx.ensure_object(dict)


def get_base(ctx):
    instance = int(ctx.parent.params['instance'])
    time = ctx.parent.params['time']
    output = ctx.parent.params['output']
    return instance, time, output


@cli.command(name='aco')
@click.pass_context
@click.option('-p', '--populationsize', required=False, type=int, default=100)
@click.option('-a', '--alpha', required=False, type=int, default=1)
@click.option('-b', '--beta', required=False, type=int, default=3)
@click.option('-r', '--rho', required=False, type=float, default=0.8)
@click.option('-l', '--localsearch', required=False, type=bool, default=False)
def ant_colony_optimization_cmd(ctx, populationsize, alpha, beta, rho, localsearch):
    click.echo('ant_colony_optimization')
    instance, time, output = get_base(ctx)
    algorithm = aco(populationsize, alpha, beta, rho, localsearch)
    algorithm.solve(instance, time)
    algorithm.report(output)


@cli.command(name='ga')
@click.pass_context
@click.option('-p', '--populationsize', required=False, type=int, default=300)
@click.option('-m', '--improvingmutation', required=False, type=bool, default=True)
@click.option('-r', '--mutationratefactor', required=False, type=int, default=1)
@click.option('-s', '--selection', required=False, type=click.Choice(['Tournament', 'SUS', 'Roulette', 'Best']), default='Tournament')
@click.option('-k', '--selectionk', required=False, type=int, default=3)
@click.option('-l', '--localsearch', required=False, type=bool, default=False)
def genetic_algorithm_cmd(ctx, populationsize, improvingmutation, mutationratefactor, selection, selectionk, localsearch):
    click.echo('genetic_algorithm')
    instance, time, output = get_base(ctx)

    algorithm = genetic_algorithm(populationsize, improvingmutation, mutationratefactor, selection, selectionk, localsearch)
    algorithm.solve(instance, time)
    algorithm.report(output)


@cli.command(name='dnla')
@click.pass_context
@click.option('-l', '--localsearch', required=False, type=bool, default=False)
def baseline_dnala_heuristic_cmd(ctx, localsearch):
    click.echo('baseline_dnala_heuristic')
    instance, time, output = get_base(ctx)
    algorithm = dnala(localsearch)
    algorithm.solve(instance, time)
    algorithm.report(output)


if __name__ == '__main__':
    cli(obj={})
