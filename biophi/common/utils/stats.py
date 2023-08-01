import hashlib
from typing import List

from flask import current_app, request
from sqlalchemy import create_engine
import os
import pandas as pd
import datetime
import time

from sqlalchemy.exc import OperationalError

from biophi.common.utils.io import AntibodyInput


def get_stats(table: str='access_log', past_days: int=7) -> pd.DataFrame:
    e = get_engine()
    from_date = (datetime.datetime.now() - datetime.timedelta(days=past_days)).date()
    from_timestamp = time.mktime(from_date.timetuple())
    try:
        stats = pd.read_sql(f'SELECT * FROM {table} WHERE timestamp >= {from_timestamp}', e, parse_dates=['timestamp'])
    except OperationalError as e:
        raise ValueError('Stats tables are not filled, try submitting a task first.', e)
    stats = stats.set_index('timestamp')
    return stats


def log_submission(antibody_inputs: List[AntibodyInput], invalid_names: List[str], duplicate_names: List[str], unrecognized_files: List[str]) -> None:
    if not current_app.config['STATS_DB_PATH']:
        return
    num_inputs_unpaired = sum((not i.heavy_protein_seq or not i.light_protein_seq) for i in antibody_inputs)
    unrecognized_extensions = set([f.split('.')[-1] for f in unrecognized_files])
    log_data(
        data=dict(
            input_hash=hashlib.md5(','.join([str(i) for i in antibody_inputs]).encode()).hexdigest(),
            num_inputs=len(antibody_inputs),
            num_inputs_unpaired=num_inputs_unpaired,
            num_inputs_heavy=sum(bool(i.heavy_protein_seq) for i in antibody_inputs),
            num_inputs_light=sum(bool(i.light_protein_seq) for i in antibody_inputs),
            num_invalid_seqs=len(invalid_names),
            num_duplicate_names=len(duplicate_names),
            unrecognized_extensions=','.join(unrecognized_extensions) if unrecognized_files else None
        ),
        table='submission_log'
    )


def log_task_result(running_seconds=None, exception=None) -> None:
    if not current_app.config['STATS_DB_PATH']:
        return
    log_data(
        data=dict(
            running_seconds=running_seconds,
            exception=str(exception) if exception else None
        ),
        table='task_log'
    )


def log_access(exception=None) -> None:
    if not current_app.config['STATS_DB_PATH']:
        return
    log_data(
        data=dict(
            method=request.method,
            path=request.path,
            exception=str(exception) if exception else None
        ),
        table='access_log'
    )


def log_data(data, table):
    start_time = time.time()
    e = get_engine()
    try:
        ip = request.headers.get('X-Forwarded-For', request.remote_addr)
        request_data = dict(
            referer=request.headers.get('Referer'),
            browser_hash=hashlib.md5(
                '{}_{}'.format(request.headers.get('User-Agent'), ip).encode()
            ).hexdigest(),
            ip=ip,
            endpoint_name=request.endpoint.split('.')[-1] if request.endpoint else None
        )
    except RuntimeError:
        request_data = {}

    df = pd.DataFrame([dict(
        timestamp=time.time(),
        **request_data,
        **data
    )]).set_index('timestamp')

    try:
        df.to_sql(table, e, if_exists='append')
    except OperationalError as e:
        print('Unable to log stats data:', e)

    duration_seconds = time.time() - start_time
    if duration_seconds > 0.5:
        print(f'Logging stats took longer than expected ({duration_seconds:.1f}s)')


_STATS_ENGINE = None


def get_engine():
    global _STATS_ENGINE
    if _STATS_ENGINE is None:
        db_path = current_app.config['STATS_DB_PATH']
        if db_path is None:
            raise ValueError('Configure STATS_DB_PATH env var to use statistics')
        _STATS_ENGINE = create_engine(
            'sqlite:///' + os.path.abspath(db_path),
            connect_args={'timeout': 5}
        )
    return _STATS_ENGINE
